import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Arrays;

public class CancerPrediction {
    private double[] thetaPreds;
    private Double[] maxDens;
    private Integer[] typeRanks;
    private double bestTheta;
    private double bestDens;
    private double bestRatio;

    public CancerPrediction(RealMatrix dens, RealVector thetas, boolean cancerPos) {
        int nTheta = dens.getColumnDimension();
        int nType = dens.getRowDimension();
        thetaPreds = new double[nType];
        maxDens = new Double[nType];
        double[] baseDens = dens.getColumn(0); // the likelihoods if the sample is normal

        for (int i=0; i<nType; i++) {
            RealVector rowDens;
            if (cancerPos) {
                rowDens = dens.getRowVector(i).getSubVector(1,nTheta - 1);
                thetaPreds[i] = thetaPred(rowDens, thetas.getSubVector(1,nTheta-1));
            }
            else {
                rowDens = dens.getRowVector(i);
                thetaPreds[i] = thetaPred(rowDens, thetas);
            }
            maxDens[i] = rowDens.getMaxValue();
        }
        //sort type by maxDens-baseDens
        typeRanks  = new Integer[nType];
        for(int i=0 ; i < nType; i++ ) typeRanks[i] = i;
        Arrays.sort(typeRanks, (i1, i2) -> { // descending
            return Double.compare(maxDens[i2]-baseDens[i2], maxDens[i1]-baseDens[i1]); // log scale
        });
        bestTheta = thetaPreds[typeRanks[0]];
        bestDens = maxDens[typeRanks[0]];
        bestRatio = bestDens-baseDens[typeRanks[0]]; // log scale
    }


    private static double thetaPred(RealVector sumLogDens, RealVector thetas) {
        int thetaIndex = sumLogDens.getMaxIndex();
        double pred;
        if (thetaIndex == -1 || // NaN for all theta dists
                sumLogDens.getMaxValue() == sumLogDens.getMinValue()) { // uniform dist
            pred = Double.NaN;
        } else {
            pred = thetas.getEntry(thetaIndex);
        }
        return pred;
    }

    public double getBestTheta() {
        return bestTheta;
    }

    public void setBestTheta(double bestTheta) {
        this.bestTheta = bestTheta;
    }

    public double[] getThetaPreds() {
        return thetaPreds;
    }

    public void setThetaPreds(double[] thetas) {
        this.thetaPreds = thetas;
    }

    public Double[] getMaxDens() {
        return maxDens;
    }

    public void setMaxDens(Double[] maxDens) {
        this.maxDens = maxDens;
    }

    public Integer[] getTypeRanks() {
        return typeRanks;
    }

    public void setTypeRanks(Integer[] typeRanks) {
        this.typeRanks = typeRanks;
    }
    public double getBestDens() {
        return bestDens;
    }

    public void setBestDens(double bestDens) {
        this.bestDens = bestDens;
    }

    public double getBestRatio() {
        return bestRatio;
    }

    public void setBestRatio(double bestRatio) {
        this.bestRatio = bestRatio;
    }
}
