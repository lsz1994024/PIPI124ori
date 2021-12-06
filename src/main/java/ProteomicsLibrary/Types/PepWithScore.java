package ProteomicsLibrary.Types;

import ProteomicsLibrary.Types.SparseBooleanVector;

public class PepWithScore implements Comparable<PepWithScore> {
    public String pepSeq;
    //    public final String[] proteins;
    public double score;
    public double count;
    public PepWithScore(String pepSeq, double score, double count) {
        this.pepSeq = pepSeq;
        this.score = score;
        this.count = count;
    }

    public int compareTo(PepWithScore other) {
        if (score < other.score) {  //decreasing order
            return 1;
        } else {
            return -1;
        }
    }


}

