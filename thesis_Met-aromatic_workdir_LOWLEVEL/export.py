"""
dsw7@sfu.ca
A script that demonstrates how data can be exported to MongoDB.
"""

from pymongo import MongoClient
from sys import path
path.append(r"utils")
from ma import MetAromatic

DEFAULT_PORT = 27017
DEFAULT_HOST = "localhost"
DB = "met_aromatic"
COL = "results"

client = MongoClient(DEFAULT_HOST, DEFAULT_PORT)
db = client[DB]
col = db[COL]

code, outgoing = '1rcy', []
for result in MetAromatic(code).met_aromatic():
    outgoing.append({
        "aro": result[0],
        "arores": result[1],
        "met": result[3],
        "norm": result[4],
        "met-theta": result[5],
        "met-phi": result[6]
    })

col.insert_many(outgoing)

for item in col.find():
    print(item)

db.drop_collection(COL)
client.close()

