/**
 *
 */
package org.theseed.proteins.kmers;

/**
 * This is a small class that is used to help find discriminating kmers.  It is used as the target of a map keyed on the
 * kmers themselves.  It tracks the role ID first associated with the kmer, the number of times the kmer was found in the
 * role, and the number of times it was found in other roles.
 *
 * @author Bruce Parrello
 *
 */
public class RoleCounter {

    // FIELDS
    /** ID of the associated role */
    private final String roleId;
    /** number of good occurrences */
    private int goodCount;
    /** number of bad occurrences */
    private int badCount;

    /**
     * Create a blank role counter.
     *
     * @param roleId	ID of the target role
     */
    public RoleCounter(String roleId) {
        this.roleId = roleId;
        this.goodCount = 0;
        this.badCount = 0;
    }

    /**
     * Denote a hit.
     *
     * @param roleHit	ID of the role that hit
     *
     * @return TRUE if the hit was good, else FALSE
     */
    public boolean count(String roleHit) {
        boolean retVal = this.roleId.contentEquals(roleHit);
        if (retVal)
            this.goodCount++;
        else
            this.badCount++;
        return retVal;
    }

    /**
     * @return TRUE if this is a good kmer
     */
    public boolean isGood() {
        return (this.badCount == 0);
    }

    /**
     * @return the ID of the relevant role
     */
    public String getRoleId() {
        return this.roleId;
    }

    /**
     * @return the number of good hits
     */
    public int getGoodCount() {
        return this.goodCount;
    }

    /**
     * @return the number of bad hits
     */
    public int getBadCount() {
        return this.badCount;
    }

}
