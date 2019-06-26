"""
dsw7@sfu.ca
The analysis script that processes the results of the Met-aromatic algorithm
"""

from pymongo import MongoClient
from networkx import Graph, connected_components
from argparse import ArgumentParser, RawTextHelpFormatter

DEFAULT_PORT = 27017
DEFAULT_HOST = "localhost"
DB = "ma"
COL = "ma"

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


def move_pairs():
    query = [
        {
            '$group': {
                '_id': {
                    'code': "$code",
                    'EC': "$ec",
                },
                'pairs': {
                    '$addToSet':  {
                        '$concat': ["$aro", "", "$arores", "|", "MET", "", "$met"]
                    }
                }
            }
        },
        {'$project': {'pairs': 1, '_id': 1}},
        {'$out': 'pairs'}
    ]

    col.aggregate(query)
    print("\n -- Wrote pair data to collection: pairs")


def get_n_bridges_from_pairs(n=3):
    if n < 3:
        exit("Incorrect bridge order. A bridge must be of n >= 3.")
    else:
        print(' -- Collected all {}-bridges: \n'.format(n))
        bridges = []
        for entry in list(db['pairs'].find()):
            pairs = []
            for pair in entry.get('pairs'):
                pairs.append(tuple(pair.split('|')))  # 'TYR123|MET123' -> ('TYR123', 'MET123')

            G = Graph()
            G.add_edges_from(pairs)

            for disconnects in list(connected_components(G)):
                if len(disconnects) == n:
                    code_and_ec = entry.get('_id')
                    code_and_ec.update({'bridge': disconnects})
                    bridges.append(code_and_ec)

        return bridges  # export back into MongoDB?


if __name__ == '__main__':
    count_entries()
    breakdowns_by_order()
    move_pairs()

    for item in get_n_bridges_from_pairs():
        print(item)


client.close()
