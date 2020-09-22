/**
 *
 */
package org.theseed.genome.compare;

import java.security.NoSuchAlgorithmException;

/**
 * This enum represents the types of comparisons that can be performed by the GenomeCompareProcessor.
 *
 * @author Bruce Parrello
 *
 */
public enum CompareType {
    FUNCTIONS, SUBSYSTEMS;

    public IGenomeMatcher create() throws NoSuchAlgorithmException {
        IGenomeMatcher retVal = null;
        switch (this) {
        case FUNCTIONS :
            retVal = new CompareGenomes();
            break;
        case SUBSYSTEMS:
            retVal = new CompareSubsystems();
            break;
        }
        return retVal;
    }

}
