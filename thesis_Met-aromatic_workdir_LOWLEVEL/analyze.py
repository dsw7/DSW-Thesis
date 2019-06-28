"""
dsw7@sfu.ca
The analysis script that processes the results of the Met-aromatic algorithm
"""

from pymongo import MongoClient, errors
from networkx import Graph, connected_components
from argparse import ArgumentParser, RawTextHelpFormatter

DEFAULT_PORT = 27017
DEFAULT_HOST = "localhost"
DB = "ma"
COL = "ma"
COLLECTION_BRIDGES = "bridges"   #TODO: need to reorganize the three primary collections
COLLECTION_PAIRS = "pairs"

msg_db = 'Choose a MongoDB database name to process. \nDefault = ma. \nUsage: $ python analyze.py --database <...>'
msg_col = 'Choose a MongoDB collection name to process. \nDefault = ma. \nUsage: $ python analyze.py --collection <...>'
msg_port = 'Choose a MongoDB port. \nDefault = 27017. \nUsage: $ python analyze.py --mongoport <port>'
msg_host = 'Choose a MongoDB host. \nDefault = localhost. \nUsage: $ python analyze.py --mongohost <host>'

parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--mongoport', help=msg_port, default=DEFAULT_PORT, type=int)
parser.add_argument('--mongohost', help=msg_host, default=DEFAULT_HOST, type=str)
parser.add_argument('--database', help=msg_db, default=DB, type=str)
parser.add_argument('--collection', help=msg_col, default=COL, type=str)

mongoport = parser.parse_args().mongoport
mongohost = parser.parse_args().mongohost
database = parser.parse_args().database
collection = parser.parse_args().collection

client = MongoClient(mongohost, mongoport)
db = client[database]
col = db[collection]

print(' -- Connected to MongoDB on:')
print(' -- Port: {}'.format(mongoport))
print(' -- Host: {}'.format(mongohost))
print(' -- Database: {}'.format(database))
print(' -- Collection: {}\n'.format(collection))
print(' -- Analyzing...\n')


def count_entries():
    count = len(list(col.distinct('code')))
    if not count:
        exit(' -- Empty database. Exiting.')
    else:
        print(' -- Number of unique entries analyzed: {}\n'.format(count))


def breakdowns_by_order():
    """
    :return: Nothing. Function performs operation pass-by style.
    """
    query = [
        {
            '$group': {
                '_id': {
                    'code': "$code",
                    'arores': "$arores",
                    'met': "$met"
                },
                'order': {
                    '$sum': 1
                }
            }
        },
        {
            '$project': {
                'order': 1
            }
        },
        {
            '$group': {
                '_id': '$order',
                'count': {
                    '$sum': 1
                },
            }
        },
        {
            '$sort': {
                '_id': 1
            }
        }
    ]

    print(' -- Breakdowns by order: ')
    cursor = col.aggregate(query)
    for entry in list(cursor):
        print(' -- Order: {} | Count: {}'.format(entry.get('_id'), entry.get('count')))


def move_pairs(name_collection='pairs'):
    """
    :param name_collection: The name of the MongoDB collection to export to.
    :return: Nothing. Function performs operation pass-by style.
    """
    query = [
        {
            '$group': {
                '_id': {
                    'code': "$code",
                    'EC': {'$substr': ["$ec", 0, 1]}  # 1.2.3.56 -> 1
                },
                'pairs': {
                    '$addToSet':  {
                        '$concat': ["$aro", "", "$arores", "|", "MET", "", "$met"]
                    }
                }
            }
        },
        {'$project': {'pairs': 1, 'EC': '$_id.EC', 'code': '$_id.code', '_id': 0}},
        {'$out': name_collection}
    ]

    col.aggregate(query)
    print("\n -- Wrote pair data to collection: pairs")


def get_bridges_from_pairs(n=2, name_collection='bridges'):
    """
    I import from MongoDB here, take advantage of NetworkX to find n-bridges using
    network theory approaches then export back to MongoDB for storage.

    :param n: 2-bridge, 3-bridge, 4-bridge, etc.
    :param name_collection: The MongoDB collection to export data to.
    :return: Nothing. Function performs operation pass-by style.
    """
    if n < 2:
        exit("Incorrect bridge order. A bridge must be of n >= 2!")
    else:
        print(' -- Collected all {}-bridges and exported to collection: bridges \n'.format(n))
        bridges = []
        for entry in list(db['pairs'].find()):
            pairs = []
            for pair in entry.get('pairs'):
                pairs.append(tuple(pair.split('|')))  # 'TYR123|MET123' -> ('TYR123', 'MET123')

            G = Graph()
            G.add_edges_from(pairs)

            for disconnects in list(connected_components(G)):
                if len(disconnects) == n + 1:
                    if ''.join(disconnects).count('MET') == 1:  # remove inverse bridges -> MET-ARO-MET
                        bridges.append(
                            {
                                'code': entry.get('code'),
                                'EC': entry.get('EC'),
                                'bridge': list(disconnects)
                            }
                        )
        try:
            db[name_collection].insert_many(bridges)
        except errors.BulkWriteError as pymongo_exception:
            print(pymongo_exception.details['writeErrors'])


if __name__ == '__main__':
    count_entries()
    breakdowns_by_order()
    move_pairs()
    get_bridges_from_pairs()

client.close()
