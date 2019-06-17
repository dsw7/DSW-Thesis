from pymongo import MongoClient

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

# count number of unique pdb entries / exit if db is empty
cursor = col.distinct('code')
count = len(list(cursor))
if not count:
    exit(' -- Empty database. Exiting.')
else:
    print(' -- Number of unique entries analyzed: {}\n'.format(count))

# count interactions by order
query_count_order = [
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
cursor = col.aggregate(query_count_order)
for entry in list(cursor):
    print(' -- Order: {} | Count: {}'.format(entry.get('_id'), entry.get('count')))

client.close()
