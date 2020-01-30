package org.theseed.proteins.kmers.anno;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import org.theseed.locations.Location;
import org.theseed.proteins.kmers.KmerReference;

/**
 * Unit test for simple App.
 */
public class AppTest extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    // SAMPLE KMERS

    private static final String KMER1 = "SNENIQNA";
    private static final String KMER2 = "SNENMQNA";
    private static final String KMER3 = "GDQM";

    /**
     * basic kmer test
     */
    public void testKmerReferences() {
        KmerReference.setKmerSize(8);
        KmerReference k1 = new KmerReference(KMER1, "contig1", 10, "+");
        assertThat(k1.getKmer(), equalTo(KMER1));
        Location loc = k1.getLoc();
        assertThat(loc.getContigId(), equalTo("contig1"));
        assertThat(loc.getLeft(), equalTo(10));
        assertThat(loc.getLength(), equalTo(24));
        assertThat(loc.getBegin(), equalTo(10));
        assertThat(loc.getDir(), equalTo('+'));
        KmerReference k1a = new KmerReference(KMER1, "contig1", 100, "-");
        assertThat(k1a.getKmer(), equalTo(KMER1));
        assertTrue(k1.equals(k1a));
        assertThat(k1.hashCode(), equalTo(k1a.hashCode()));
        loc = k1a.getLoc();
        assertThat(loc.getContigId(), equalTo("contig1"));
        assertThat(loc.getLeft(), equalTo(100));
        assertThat(loc.getLength(), equalTo(24));
        assertThat(loc.getEnd(), equalTo(100));
        assertThat(loc.getDir(), equalTo('-'));
        KmerReference k2 = new KmerReference(KMER2, "contig1", 100, "-");
        assertThat(k2.getKmer(), equalTo(KMER2));
        assertFalse(k1.equals(k2));
        Location loc2 = k2.getLoc();
        assertThat(loc2, equalTo(loc));
        KmerReference.setKmerSize(4);
        KmerReference k3 = new KmerReference(KMER3, "contig2", 10, "+");
        assertThat(k3.getKmer(), equalTo(KMER3));
        loc = k3.getLoc();
        assertThat(loc.getContigId(), equalTo("contig2"));
        assertThat(loc.getLeft(), equalTo(10));
        assertThat(loc.getBegin(), equalTo(10));
        assertThat(loc.getLength(), equalTo(12));
    }


}
