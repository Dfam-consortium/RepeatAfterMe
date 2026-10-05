#
# RepeatAfterMe
#   Top-level build for the Rust implementation. The C sources that produced
#   releases through V0.0.7 live under c/ and build separately; see `make c`.
#
# Single source of truth for the version is [workspace.package] in Cargo.toml.
VERSION := $(shell sed -n 's/^version *= *"\(.*\)"/\1/p' Cargo.toml | head -1)
INSTDIR = /usr/local/RepeatAfterMe-$(VERSION)
CARGO   = cargo

RUST_SRC = $(wildcard ram-core/src/*.rs ram-cli/src/*.rs) \
           Cargo.toml ram-core/Cargo.toml ram-cli/Cargo.toml

.PHONY: all install c test clean

all: RAMExtend

# Cargo names the binary ram-extend. RepeatModeler probes $(INSTDIR)/RAMExtend
# and dfam-tetools runs RAMExtend -version from the top of an unpacked release,
# so a bare `make` leaves that name here as well. Both names are the same
# binary.
RAMExtend: target/release/ram-extend
	cp -f $< $@

target/release/ram-extend: $(RUST_SRC)
	$(CARGO) build --release

install: all
	mkdir -p $(INSTDIR)
	cp target/release/ram-extend $(INSTDIR)/ram-extend
	ln -sf ram-extend $(INSTDIR)/RAMExtend
	cp README.md $(INSTDIR)
	mkdir -p $(INSTDIR)/util
	install -m 755 util/extend-stk.pl $(INSTDIR)/util/extend-stk.pl

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
