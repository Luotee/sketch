hot:
	gcc -o hotitems hotitems.c prng.c cgt.c lossycount.c massdal.c  frequent.c ccfc.c countmin.c -lm -Wall 
stable: 
	gcc -o teststab teststab.c prng.c massdal.c stable.c ams.c ccfc.c fm.c -lm -Wall
change: change.c changewrapper.c countmin.c
	gcc -o change changewrapper.c prng.c massdal.c change.c countmin.c -lm -Wall -O3
