package ProteomicsLibrary.Types;

import ProteomicsLibrary.Types.SparseBooleanVector;

import java.util.Set;

public class PepWithScore implements Comparable<PepWithScore> {
    public String pepSeq;
    //    public final String[] proteins;
    public double score;
    public double count;
    public boolean isDecoy = false;
    public Set<String> proteins;
    public boolean hasPTM;
    public PepWithScore(String pepSeq, double score, double count, boolean isDecoy, Set<String> proteins, boolean hasPTM) {
        this.pepSeq = pepSeq;
        this.score = score;
        this.count = count;
        this.isDecoy = isDecoy;
        this.proteins = proteins;
        this.hasPTM = hasPTM;
    }

    public int compareTo(PepWithScore other) {
        if (score < other.score) {  //decreasing order
            return 1;
        } else {
            return -1;
        }
    }

}

