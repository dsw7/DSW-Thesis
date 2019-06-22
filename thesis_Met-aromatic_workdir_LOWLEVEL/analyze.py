from pymongo import MongoClient
from networkx import Graph, connected_components

DEFAULT_PORT = 27017
DEFAULT_HOST = "localhost"
DB = "ma"
COL = "ma"

print(' -- Connected to MongoDB on:')
print(' -- Port: {}'.format(DEFAULT_PORT))
print(' -- Host: {}'.format(DEFAULT_HOST))
print(' -- Database: {}'.format(DB))
print(' -- Collection: {}\n'.format(COL))
print(' -- Analyzing...\n')

client = MongoClient(DEFAULT_HOST, DEFAULT_PORT)
db = client[DB]
col = db[COL]


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
        print(' -- Collected all {}-bridges: '.format(n))
        bridges = []
        for entry in list(db['pairs'].find()):
            pairs = []
            for pair in entry.get('pairs'):
                pairs.append(tuple(pair.split('|')))  # 'TYR123|MET123' -> ('TYR123', 'MET123')

            G = Graph()
            G.add_edges_from(pairs)

            for disconnects in list(connected_components(G)):
                if len(disconnects) == n:
                    bridges.append((entry.get('_id'), disconnects))

        return bridges


if __name__ == '__main__':
    count_entries()
    breakdowns_by_order()
    move_pairs()

    for item in get_n_bridges_from_pairs():
        print(item)


client.close()
