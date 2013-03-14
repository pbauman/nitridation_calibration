
default: all
	@true # override later % rule

dirlist := $(wildcard */)

clean:
	@echo make clean
	@for dir in $(dirlist); do \
	  if [ -f "$${dir}"/Makefile ]; then \
	    $(MAKE) -C "$${dir}" clean || touch "$${dir}"/clean_failed; \
	  fi; \
	done

%:
	@echo make $@
	@for dir in $(dirlist); do \
	  if [ -f "$${dir}"/bootstrap ]; then \
            echo "Running bootstrap in $${dir}"; \
	    (cd "$${dir}" && ./bootstrap || touch bootstrap_failed); \
	  fi; \
	  if [ -f "$${dir}"/configure ]; then \
            echo "Running configure in $${dir}"; \
	    (cd "$${dir}" && ./configure || touch bootstrap_failed); \
	  fi; \
	  if [ -f "$${dir}"/Makefile ]; then \
	    $(MAKE) -C "$${dir}" $@ || touch "$${dir}"/make_failed; \
	  fi; \
	done
