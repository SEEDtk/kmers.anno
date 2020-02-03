/**
 *
 */
package org.theseed.proteins.kmers;

import java.util.Collection;
import java.util.Map;

import org.theseed.genome.Genome;
import org.theseed.locations.Location;

/**
 * This is the base class for extracting the contig kmer map from a genome.  Its
 * purpose is to take a genome as input and return a map of kmers to locations.
 *
 * @author Bruce Parrello
 *
 */
public abstract class KmerFactory {

    /** types of kmer factors */
    public static enum Type {
        /** only return unique kmers */
        STRICT,
        /** return all kmers */
        AGGRESSIVE
        ;
    }

    /**
     * @return a kmer factory of the proposed type
     *
     * @param type	type of kmer factory to create
     */
    public static KmerFactory create(KmerFactory.Type type) {
        KmerFactory retVal = null;
        switch (type) {
        case STRICT :
            retVal = new KmerFactory.Strict();
            break;
        case AGGRESSIVE :
            retVal = new KmerFactory.Aggressive();
            break;
        default :
            throw new IllegalArgumentException("Unsupported contig kmer algorithm.");
        }
        return retVal;
    }

    /**
     * @return a map of kmers to locations in the specified genome
     *
     * @param genome	genome of interest
     */
    public abstract Map<String, Collection<Location>> findKmers(Genome genome);


    /**
     * The STRICT method only returns unique kmers.
     */
    public static class Strict extends KmerFactory {

        @Override
        public Map<String, Collection<Location>> findKmers(Genome genome) {
            Map<String, Collection<Location>> retVal = KmerReference.getContigKmers(genome);
            retVal.entrySet().removeIf(x -> (x.getValue().size() > 1));
            return retVal;
        }

    }

    /**
     * The AGGRESSIVE method returns all kmers.
     */
    public static class Aggressive extends KmerFactory {

        @Override
        public Map<String, Collection<Location>> findKmers(Genome genome) {
            return KmerReference.getContigKmers(genome);
        }

    }



}
