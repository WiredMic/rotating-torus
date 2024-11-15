##
# Rotaing Donut
#
# @file
# @version 0.1

all: clean build run

clean:
	rm -f donut

build:
	# When building include all files
	# Include math
	gcc donut.c vector.c bivector.c mvec.c rotor.c -o donut -lm

run:
	./donut

# end
