package ProteomicsLibrary.Types;

//import it.unimi.dsi.fastutil.Hash;
//import sun.print.SunPageSelection;

import java.util.*;

public class SparseBooleanVector {
    // this is peptide coded vector

    public HashMap<Integer, Double> sparseVector = new HashMap<>();

    public double get(int i) {
        if (sparseVector.containsKey(i)) {
            return sparseVector.get(i);
        } else {
            return 0;
        }
    }

    public Set<Integer> getKeys() {
        return sparseVector.keySet();
    }

    public SparseBooleanVector(HashMap<Integer, Double> tag_idx) {
        sparseVector = tag_idx;
//        this.sparseVector = new HashSet<>(sparseVector);
    }

    public SparseBooleanVector() {
    }

    public void add(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            if (sparseVector.containsKey(i)) {
                sparseVector.put(i, sparseVector.get(i) + v);
            } else {
                sparseVector.put(i, v);
            }
        }
    }

    public void put(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            sparseVector.put(i, v);
        }
    }

    public double norm2square() {
        double output = 0;
        for (double v : sparseVector.values()) {
            output += v * v;
        }
        return output;
    }

    public boolean compare(SparseBooleanVector other){
        boolean flag = true;
        Map<Integer, Double> otherVector = other.sparseVector;

        if(sparseVector.size() != otherVector.size()){
            return false;
        }

        for (int i : sparseVector.keySet()) {
            if (!(sparseVector.get(i).equals(otherVector.get(i)))) {
                flag = false;
            }
        }

        for (int i : otherVector.keySet()) {
            if (!(otherVector.get(i).equals(sparseVector.get(i)))){
                flag = false;
            }
        }
        return flag;
    }

//    public double dot(SparseBooleanVector other) {
//        double output = 0;
//        for (int i : sparseVector.keySet()) {
//            output += other.get(i);
//        }
//        return output;
//    }

    public double dot(SparseBooleanVector other) {
        double output = 0;
        Map<Integer, Double> otherVector = other.sparseVector;
        Set<Integer> intersectedKeys = new HashSet<>(sparseVector.keySet());
        intersectedKeys.retainAll(otherVector.keySet());
        for (int i : intersectedKeys) {
            output += sparseVector.get(i) * otherVector.get(i);
        }
        return output;
    }

//    public List sort(SparseBooleanVector vector){
//        List<Map.Entry<Integer, Double>> list_Data = new ArrayList<Map.Entry<Integer, Double>>(vector.sparseVector.entrySet());
//        Collections.sort(list_Data, new Comparator<Map.Entry<Integer, Double>>()
//        {
//            public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2)
//            {
//                if ((o2.getValue() - o1.getValue())>0)
//                    return 1;
//                else if((o2.getValue() - o1.getValue())==0)
//                    return 0;
//                else
//                    return -1;
//            }
//        });
//        return list_Data();
//    }

    public SparseBooleanVector normalization(SparseBooleanVector code){
        Set<Integer> newset = code.getKeys();
        for (int i : newset){
            code.put(i, 1.0);
        }
        return code;
    }

    public double sqrt_cos(SparseBooleanVector other){
        double sum1 = 0, sum2 = 0;
        Map<Integer, Double> otherVector = other.sparseVector;
        for(double temp : sparseVector.values()){
            sum1 += temp;
        }
        for(double temp : otherVector.values()){
            sum2 += temp;
        }

        return sum1 * sum2;
    }

    public double intersection(SparseBooleanVector other){
        double output = 0;
        Map<Integer, Double> otherVector = other.sparseVector;
        Set<Integer> intersectedKeys = new HashSet<>(sparseVector.keySet());
        intersectedKeys.retainAll(otherVector.keySet());
        int l1 = otherVector.size();
        int l2 = sparseVector.size();

        return (l1 > l2) ? (intersectedKeys.size() / l2) : (intersectedKeys.size() / l1);
    }

    public double pearson(SparseBooleanVector other) {
        double tempsum1 = 0;
        double tempno1 = 0, tempno2 = 0, tempde1 = 0, tempde2 = 0;
        for(double temp1 : sparseVector.values()){
            tempsum1 += temp1;
        }
        double tempmean1 = tempsum1 / sparseVector.size();

        Map<Integer, Double> otherVector = other.sparseVector;
        double tempsum2 = 0;
        for(double temp2 : otherVector.values()){
            tempsum2 += temp2;
        }
        double tempmean2 = tempsum2 / otherVector.size();

        for(double temp1 : sparseVector.values()){
            tempno1 += temp1 - tempmean1;
        }
        for(double temp2 : otherVector.values()){
            tempno2 += temp2 - tempmean2;
        }

        double nominator = tempno1 * tempno2;

        for (double temp1 : sparseVector.values()){
            tempde1 += Math.pow((tempmean1 - temp1),2);
        }

        for (double temp2 : otherVector.values()){
            tempde2 += Math.pow((tempmean2 - temp2),2);
        }

        double denominator = tempde1 * tempde2;

        return nominator / denominator;
    }

//    public double euclidean(SparseBooleanVector other) {
//        double output = 0;
//        Map<Integer, Double> otherVector = other.sparseVector;
//        Set<Integer> unionKeys = new HashSet<>();
//        unionKeys.addAll(sparseVector.keySet());
//        unionKeys.addAll(otherVector.keySet());
//        for (int i : unionKeys) {
//            if(sparseVector.containsKey(i) && otherVector.containsKey(i)){
//                output += Math.pow((sparseVector.get(i) - otherVector.get(i)), 2);
//            }
//            else if (sparseVector.containsKey(i) && !otherVector.containsKey(i)){
//                output += Math.pow(sparseVector.get(i), 2);
//            }
//            else if (!sparseVector.containsKey(i) && otherVector.containsKey(i)){
//                output += Math.pow(otherVector.get(i), 2);
//            }
//        }
//        return output / unionKeys.size();
//    }

    public double getMaxValue() {
        List<Double> intensityList = new ArrayList<>(sparseVector.values());
        intensityList.sort(Collections.reverseOrder());
        return intensityList.get(0);
    }


//    public double fastDot(SparseVector other) {
//        double output = 0;
//        Map<Integer, Double> otherVector = other.getVectorMap();
//        sparseVector.retainAll(otherVector.keySet());
//        for (int i : sparseVector.keySet()) {
//            output += otherVector.get(i);
//        }
//        return output;
//    }

//    public double dot(SparseBooleanVector other) {
//        Set<Integer> intersectedKeys = new HashSet<>(sparseVector);
//        intersectedKeys.retainAll(other.sparseVector);
//        return intersectedKeys.size();
//    }

//    public SparseBooleanVector deepCopy() {
//        return new SparseBooleanVector(this.sparseVector);
//    }

//    public boolean isZero(int idx) {
//        return !sparseVector.contains(idx);
//    }

//    public void delete(int idx) {
//        if (sparseVector.contains(idx)) {
//            sparseVector.remove(idx);
//        }
//    }

    public int getNonZeroNum() {
        return sparseVector.size();
    }

//    public Integer[] getNonZeroIdxes() {
//        return sparseVector.toArray(new Integer[0]);
//    }

    public String toString() {
        StringBuilder sb = new StringBuilder(sparseVector.size() * 6);
        for (int idx : sparseVector.keySet()) {
            sb.append(idx);
            sb.append(";");
        }
        return sb.toString();
    }
}
