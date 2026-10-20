#
# RepeatAfterMe
#   Top-level build for the Rust implementation. The C sources that produced
#   releases through V0.0.7 live under c/ and build separately; see `make c`.
#
# Layout produced by `make install` (flat, because RepeatModeler looks for
# $REPEATAFTERME_DIR/RAMExtend directly):
#   $(PREFIX)/$(NAME)-$(VERSION)/ram-extend           the compiled binary
#   $(PREFIX)/$(NAME)-$(VERSION)/RAMExtend            symlink to ram-extend
#   $(PREFIX)/$(NAME)-$(VERSION)/util/extend-stk.pl   finds ../ram-extend
#   $(PREFIX)/$(NAME)-$(VERSION)/README.md, LICENSE   documentation (see DOCS)
#
# Examples:
#   make                          # build, and copy RAMExtend to the top of the tree
#   make install PREFIX=/opt      # -> /opt/RepeatAfterMe-<version>/...
#   make dist PREFIX=/tmp/stage   # build a versioned tarball
#   make version                  # print the workspace version from Cargo.toml

PREFIX   ?= /usr/local
DESTDIR  ?=
NAME     ?= RepeatAfterMe
CARGO    ?= cargo
INSTALL  ?= install
MANIFEST  = Cargo.toml
# Documentation shipped in a binary release. Listed by name rather than globbed
# as *.md, so that a maintainer note at the top level (TODO.md, CHANGELOG.md)
# does not silently end up in a release tarball.
DOCS      = README.md LICENSE

# Single source of truth for the version is [workspace.package] in Cargo.toml.
# The sed range confines the match to that table.
VERSION := $(shell sed -n '/^\[workspace\.package\]/,/^\[/s/^version *= *"\([^"]*\)".*/\1/p' $(MANIFEST) | head -1)
$(if $(VERSION),,$(error could not read version from $(MANIFEST)))

# Versioned package root (DESTDIR-aware for staged packaging). INSTDIR may be
# set directly to install somewhere that does not follow the naming.
INSTDIR  ?= $(DESTDIR)$(PREFIX)/$(NAME)-$(VERSION)

RUST_SRC = $(wildcard ram-core/src/*.rs ram-cli/src/*.rs) \
           Cargo.toml ram-core/Cargo.toml ram-cli/Cargo.toml

.PHONY: all version install uninstall dist c test clean

all: RAMExtend

version:
	@echo $(VERSION)

# Cargo names the binary ram-extend. RepeatModeler probes $(INSTDIR)/RAMExtend
# and dfam-tetools runs RAMExtend -version from the top of an unpacked release,
# so a bare `make` leaves that name here as well. Both names are the same
# binary.
RAMExtend: target/release/ram-extend
	cp -f $< $@

target/release/ram-extend: $(RUST_SRC)
	$(CARGO) build --release

install: all
	$(INSTALL) -d $(INSTDIR)
	$(INSTALL) -m 0755 target/release/ram-extend $(INSTDIR)/ram-extend
	ln -sf ram-extend $(INSTDIR)/RAMExtend
	@for d in $(DOCS); do \
	  $(INSTALL) -m 0644 $$d $(INSTDIR)/ && echo "  doc  <- $$d"; \
	done
	$(INSTALL) -d $(INSTDIR)/util
	$(INSTALL) -m 0755 util/extend-stk.pl $(INSTDIR)/util/extend-stk.pl
	@echo "Installed $(NAME) $(VERSION) -> $(INSTDIR)"

uninstall:
	rm -rf $(INSTDIR)

# Versioned tarball of the installed hierarchy: <name>-<version>.tar.gz
dist: install
	tar -C $(DESTDIR)$(PREFIX) -czf $(NAME)-$(VERSION).tar.gz $(NAME)-$(VERSION)
	@echo "Wrote $(NAME)-$(VERSION).tar.gz"

# The frozen C tree. Its Makefile is unchanged from V0.0.7.
c:
	$(MAKE) -C c

test:
	$(CARGO) test
	harness/diff-c-rust.sh

clean:
	$(CARGO) clean
	-rm -f RAMExtend
	$(MAKE) -C c clean
