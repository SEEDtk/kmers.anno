/**
 *
 */
package org.theseed.proteins.kmers;

import org.theseed.genome.Contig;

/**
 * This class computes the position of a kmer in a specified strand.  Given a position in the strand's DNA
 * string, it returns the real left position of the kmer.  For the Plus strand this is the position itself.
 * For the minus strand we are looking at a reversed sequence, so we need to subtract from the length of
 * the contig.
 *
 * @author Bruce Parrello
 *
 */
public abstract class KmerPosition {

    // FIELDS
    /* length of the contig */
    protected int contigLen;
    /* length of the kmers */
    protected int kmerLen;

    /**
     * Construct a blank Kmer Position object.
     *
     * @param contig	the target contig
     */
    public KmerPosition(Contig contig) {
        this.contigLen = contig.length();
        this.kmerLen = KmerReference.getKmerSize() * 3;
    }

    /**
     * Compute the left edge of the kmer found at the specified position in the current contig sequence.
     *
     * @param pos	the offset into the translated protein string where the kmer was found
     * @param frame current frame offset in the contig
     *
     * @return the contig position for the left edge of the kmer
     */
    public abstract int calcLeft(int pos, int frame);

    /**
     * @return the current strand
     */
    public abstract String getDir();

    /**
     * Plus strand subclass.
     */
    public static class Plus extends KmerPosition {

        public Plus(Contig contig) {
            super(contig);
        }

        @Override
        public int calcLeft(int pos, int frame) {
            return pos * 3 + frame;
        }

        @Override
        public String getDir() {
            return "+";
        }
    }

    /**
     * Minus strand subclass.
     */
    public static class Minus extends KmerPosition {

        /** base number used to convert the position to a left edge */
        int base;

        public Minus(Contig contig) {
            super(contig);
            this.base = contigLen - kmerLen + 2;
        }

        @Override
        public int calcLeft(int pos, int frame) {
            return this.base - (pos * 3 + frame);
        }

        @Override
        public String getDir() {
            return "-";
        }

    }
}
