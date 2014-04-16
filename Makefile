LIB_DIR := lib
BIN_DIR := bin
APP	:= $(BIN_DIR)/mwiggle
LD_FLAGS	:= -Llib -lmw 

all:	$(APP)

test: $(APP)
	cd tests; $(MAKE)

$(APP): check
	cd src; $(MAKE)


check: $(LIB_DIR) $(BIN_DIR)

$(LIB_DIR):
	@mkdir $@

$(BIN_DIR):
	@mkdir $@

clean: 
	@rm -rf $(BIN_DIR) $(LIB_DIR)
	cd src; $(MAKE) clean
	cd tests; $(MAKE) clean


