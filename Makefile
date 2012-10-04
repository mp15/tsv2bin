tsv2bin: tsv2bin.c
	gcc -g -std=c99 -Wall -o tsv2bin tsv2bin.c
test: tsv2bin
	./tsv2bin -m test_data/CD_Rep1441448.map -i test_data/CD_Rep1441448.csv -o test_data/CD_Rep1441448.GTypeX
clean:
	rm tsv2bin tsv2bin.dSYM

