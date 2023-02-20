# Create the variables needed.
SCAFF=17

rm -f tmp.db && awk 'BEGIN {printf("create table T(Position,Rate); \
	BEGIN TRANSACTION;\n");}{printf("INSERT INTO T(Position,Rate) VALUES(\"%s\",%s);\n",$1,$2);} \
	END {printf("COMMIT; SELECT (Position/1000000)*1000000 as G,((Position/1000000)+1)*1000000,AVG(Rate) FROM T GROUP BY G;");}' \
	scaffold-$SCAFF-positions-own-2.txt | sqlite3 -separator $'\t' tmp.db > scaffold-$SCAFF-1Mb-2.txt && rm -f tmp.db
