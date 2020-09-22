/**
 *
 */
package org.theseed.genome.compare;

import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;
import java.util.Set;
import java.util.stream.Collectors;

import org.theseed.genome.Genome;
import org.theseed.genome.SubsystemRow;

/**
 * This class compares the subsystems in sequence-identical genomes.  A subsystem in the new genome is considered good if
 * it exists in the old genome and bad otherwise.  The new genome is considered correct, and we are downgrading the old genome
 * by the number of false negatives.
 *
 * @author Bruce Parrello
 *
 */
public class CompareSubsystems extends MatchGenomes implements IGenomeMatcher {

    /**
     * Construct a subsystem-comparison engine.
     *
     * @throws NoSuchAlgorithmException
     */
    public CompareSubsystems() throws NoSuchAlgorithmException {
        super();
    }

    // FIELDS
    /** number of good subsystems */
    private int good;
    /** number of bad subsystems */
    private int bad;

    @Override
    public boolean compare(Genome newGenome, Genome oldGenome) throws UnsupportedEncodingException {
        this.good = 0;
        this.bad = 0;
        // Get all the subsystems in the old genome.
        Set<String> oldSubs = oldGenome.getSubsystems().stream().map(x -> x.getName()).collect(Collectors.toSet());
        // Count the subsystems in the new genome.
        for (SubsystemRow sub : newGenome.getSubsystems()) {
            String subName = sub.getName();
            if (oldSubs.contains(subName))
                this.good++;
            else
                this.bad++;
        }
        // There is no problem with genome mismatch here, so we always return TRUE.
        return true;
    }

    @Override
    public int getGood() {
        return this.good;
    }

    @Override
    public int getBad() {
        return this.bad;
    }

    @Override
    public double percent() {
        double retVal = 0.0;
        if (this.good > 0)
            retVal = (good * 100.0) / (good + bad);
        return retVal;
    }
    // FIELDS
    // TODO data members for CompareSubsystems

    // TODO constructors and methods for CompareSubsystems
}
