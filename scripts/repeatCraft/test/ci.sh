#!/bin/bash
#strict mode
python repeatcraft.py -r example/example_input.gff -u example/example_input.out -c example/repeatcraft_strict.cfg -o example_strict -m strict

# loose mode
python repeatcraft.py -r example/example_input.gff -u example/example_input.out -c example/repeatcraft.cfg -o example_loose -m loose


