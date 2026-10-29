release:
	make -f release.makefile

debug:
	make -f debug.makefile

run: release
	./aiar.out

.PHONY: clean
clean:
	rm -f *.o *.mod *.MOD *genmod* *~ aiar.out aiar_dbg.out
