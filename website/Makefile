
FILES=documentation downloads examples htslib index license links man_latest man_0112b man_0112a perl_examples perl_module specs

install:
	    @for i in $(FILES); do echo "php src/index.php $$i > $$i.html"; php src/index.php $$i > $$i.html; done;
