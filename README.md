# LongAGE: defining breakpoints of genomic structural variants through optimal and memory efficient alignments of long reads

LongAGE is a memory-efficient implementation of AGE (https://github.com/abyzovlab/AGE). Its memory footprint is less than hundreds of megabytes, while it is at most four times slower than AGE  in  terms of running time.  

The tool refines the resolution and standardization of structural variant (SV) breakpoints in highly repetitive regions at a single base pair. It is capable of refining read alignment once a read has been heuristically mapped to a particular genomic location that is expected to contain an SV.

## Compilation
```
git clone https://github.com/Coaxecva/LongAGE.git
cd LongAGE
make
```

## Running
```
./long_age_align file1.fa file2.fa
```

## Help
```
./long_age_align
```

## Examples
```
# INDEL mode (default)
./long_age_align part_chr19.fa seq_first.fa 
./long_age_align -indel part_chr19.fa seq_first.fa 
# TDUP mode
./long_age_align -tdup part_chr19.fa seq_second.fa
```

## Bugs

Feel free to report any bugs that you find on GitHub:
https://github.com/Coaxecva/LongAGE/issues

## License

Released under MIT licence.
