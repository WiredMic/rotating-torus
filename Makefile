##
# Rotaing Donut
#
# @file
# @version 0.1

all: clean build run

clean:
	rm -f donut.bin

build:
	# When building include all files
	# Include math
	gcc donut.c vector.c bivector.c mvec.c rotor.c -o donut.bin -lm -Wall -Wextra -Werror

run:
	./donut.bin

# end
