SUBDIRS = src

EXTRA_DIST = LICENSE README.md examples

# Create a '.tarball-version' file containing the version string
# when creating a distribution tarball
# (not needed when using git repository)
dist-hook:
	@echo $(VERSION) > $(distdir)/.tarball-version
