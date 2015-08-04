
BIN = fill-aa fill-an-ac fill-fs fill-ref-md5 vcf-annotate vcf-compare vcf-concat \
    vcf-consensus vcf-contrast vcf-convert vcf-fix-ploidy vcf-indel-stats vcf-isec vcf-merge vcf-phased-join \
    vcf-query vcf-shuffle-cols vcf-sort vcf-stats vcf-subset vcf-to-tab vcf-tstv vcf-validator
MOD = FaSlice.pm Vcf.pm VcfStats.pm

.PHONY: install
install:
	for i in $(BIN); do install -Dm755 $(CURDIR)/$$i $(BINDIR)/$$i; done
	for i in $(MOD); do install -Dm644 $(CURDIR)/$$i $(MODDIR)/$$i; done

.PHONY: clean
clean:
	for i in $(BIN); do rm -f $(BINDIR)/$$i; done
	for i in $(MOD); do rm -f $(MODDIR)/$$i; done
