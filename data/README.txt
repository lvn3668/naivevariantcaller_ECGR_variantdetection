Readme

Naivevariantcaller takes as input reference and reads.
It indexes the reference over 1-3 bp as we are looking for ECGR point mutations
It finds all supporting reads moving along the reference genome starting with the lowest base pair that has been sequenced.

For that bunch of supporting reads, it calculates at each position, the variant, the dimer and the upstream/downstream trimer in keeping with problem definition.

It then checks against reference genome at the same position and classifies it as a point mutation if the number of supporting reads is greater than 5.

It writes out to csv position along reference genome, reference, alt, and counts of supporting reads.

It assumes no indels and no base calling errors.

Noise is defined as mutations supported by fewer than 5 reads.

All mutations of 1-3 bp length are classified as point mutations.
For 3 bp mutations, the check constitutes comparing trimers and finding no overlapping bases.
For 2 bp mutations, the check constitutes comparing dimers and finding no overlapping bases
The longest point mutation is written out to file, thus a unique trimer is not written out as 3 point mutations




