# Makefile for generating R packages.
# 2011 Andrew Redd, 2024 Sodbo Sharapov
#
# Assumes Makefile is in a folder where package contents are in a subfolder pkg.
# Roxygen uses the roxygen2 package, and will run automatically on check and all.

PKG_VERSION=$(shell grep -i ^version DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

build:
	R CMD build ./

all: roxygen build check install clean

check: $(PKG_NAME)_$(PKG_VERSION).tar.gz roxygen
	R CMD check $(PKG_NAME)_$(PKG_VERSION).tar.gz --no-manual

install: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD INSTALL $(PKG_NAME)_$(PKG_VERSION).tar.gz
	
roxygen:
	Rscript -e "library(roxygen2);roxygenize('./')"

clean:
	-rm -f $(PKG_NAME)_*.tar.gz
	-rm -r -f $(PKG_NAME).Rcheck
