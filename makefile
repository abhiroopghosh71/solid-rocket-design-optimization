LIB_DIR = librocket

default: approcket

approcket: setup_rocket.py approcket.pyx $(LIB_DIR)/libapprocket.a
	python3 setup_rocket.py build_ext --inplace && rm -f approcket.c && rm -Rf build

$(LIB_DIR)/libapprocket.a:
	make -C $(LIB_DIR) libapprocket.a

clean:
	rm -f *.so
	rm -f $(LIB_DIR)/*.so $(LIB_DIR)/*.a $(LIB_DIR)/*.o $(LIB_DIR)/*.out
