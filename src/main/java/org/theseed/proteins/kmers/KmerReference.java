/**
 *
 */
package org.theseed.proteins.kmers;

import org.theseed.locations.Location;

/**
 * This object represents a protein kmer at a specific location.  It is keyed on the kmer
 * text, so two identical kmers at different locations will compare equal.
 *
 * @author Bruce Parrello
 *
 */
public class KmerReference {

    /** protein kmer size */
    private static int K = 8;

    /** offset from start of kmer to end */
    private static int K3 = 23;

    // FIELDS

    /** kmer sequence */
    String kmer;
    /** location of the kmer */
    Location loc;

    /**
     * Change the kmer size.
     *
     * @param newSize	proposed new kmer size
     */
    public static void setKmerSize(int newSize) {
        K = newSize;
        K3 = K * 3 - 1;
    }

    /**
     * @return the kmer size
     */
    public static int getKmerSize() {
        return K;
    }

    /**
     * Create a reference to a specific kmer at a specific location.
     *
     * @param kmer		protein kmer sequence
     * @param seqId		ID of the sequence containing the kmer
     * @param begin		start position of the kmer (1-based)
     * @param dir		strand of the kmer (+ or -)
     */
    public KmerReference(String kmer, String seqId, int begin, String dir) {
        // Create the location.
        this.loc = Location.create(seqId, dir, begin, begin + K3);
        // Store the kmer string.
        this.kmer = kmer;
    }

    /**
     * @return the hash code for this object; note we only hash the kmer
     */
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((kmer == null) ? 0 : kmer.hashCode());
        return result;
    }

    /**
     * @return TRUE if this kmer reference is for the same kmer as another reference
     */
    @Override
    public boolean equals(Object obj) {
        boolean retVal;
        if (this == obj)
            retVal = true;
        else if (!(obj instanceof KmerReference))
            retVal = false;
        else {
            KmerReference other = (KmerReference) obj;
            if (kmer == null)
                retVal = (other.kmer == null);
            else
                retVal = kmer.equals(other.kmer);
        }
        return retVal;
    }

    /**
     * @return the kmer sequence
     */
    public String getKmer() {
        return kmer;
    }

    /**
     * @return the location of the kmer
     */
    public Location getLoc() {
        return loc;
    }


}
