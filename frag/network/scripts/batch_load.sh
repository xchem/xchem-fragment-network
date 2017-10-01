bin/neo4j-import --into /var/lib/neo4j/data/databases/graph.db --id-type string \
                     --nodes:Node customers.csv --nodes products.csv  \
                     --nodes orders_header.csv,orders1.csv,orders2.csv \
                 --relationships:CONTAINS order_details.csv \
                 --relationships:ORDERED customer_orders_header.csv,orders1.csv,orders2.csv


USING PERIODIC COMMIT LOAD CSV FROM 'file:///jm7b00809_si_002.txt ' AS line FIELDTERMINATOR
' ' MERGE (:F2 { smiles: line[1], hac: toInt(line[2]), chac: toInt(line[3]), osmiles: line[4]});
USING PERIODIC COMMIT LOAD CSV FROM 'file:///jm7b00809_si_003.txt ' AS line FIELDTERMINATOR
' ' MATCH (n1:F2 { smiles: line[1]}), (n2:F2 { smiles: line[2]}) MERGE (n1)-[:F2EDGE{label:line[3]}]-
>(n2);
USING PERIODIC COMMIT LOAD CSV FROM 'file:///jm7b00809_si_004.txt ' AS line FIELDTERMINATOR
' ' MATCH (n:F2 { smiles: line[1]} ) set n:MOL, n:EM, n.EM=toInt(line[3]);