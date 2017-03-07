import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Stream;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.FastMath;

/**
 * Created by Shuli Kang on 3/3/15.
 */
public class MethyModel {
	private RealVector alpha;
	private RealVector beta;
	private RealVector mean;
	private RealVector sd;
	private RealVector naRatio;
	private int nSample;
	private int nFeature;

	public MethyModel(RealMatrix m) {
		nFeature = m.getColumnDimension();
		nSample = m.getRowDimension();
		alpha = new ArrayRealVector(nFeature);
	    beta = new ArrayRealVector(nFeature);
		mean = new ArrayRealVector(nFeature);
		sd = new ArrayRealVector(nFeature);
	    naRatio = new ArrayRealVector(nFeature);
		for (int i = 0; i < nFeature; i++) {
			double [] values = m.getColumn(i);
			double sum=0;
		    List<Double> nonNa = new ArrayList<>();
			for (int j = 0; j<nSample; j++) {
				if (!Double.isNaN(values[j])) {
					sum += values[j];
					nonNa.add(values[j]);
				}
			}
			int nonNaCount = nonNa.size();

			double[] para = {Double.NaN, Double.NaN, Double.NaN, Double.NaN};
			if (nonNaCount >=2) {
				para = estBetaDist(Stream.of(nonNa.toArray(new Double[nonNaCount])).
						mapToDouble(Double::doubleValue).toArray());
			}

			alpha.setEntry(i, para[0]);
			beta.setEntry(i, para[1]);
			mean.setEntry(i,para[2]);
			sd.setEntry(i,para[3]);
			naRatio.setEntry(i, (double)(nSample-nonNaCount)/nSample);
		}
	}
	
	public MethyModel() {
	}


	public MethyModel selectFeature(boolean[] selectedFeatures) {
		MethyModel newModel = new MethyModel();
		int maxFeature = selectedFeatures.length;
		double[] newAlpha = new double[maxFeature];
		double[] newBeta = new double[maxFeature];
		double[] newMean = new double[maxFeature];
		double[] newSd = new double[maxFeature];
		double[] newNaRatio = new double[maxFeature];
		int j = 0; // num of selected features
		for(int i=0; i<maxFeature; i++) {
			if (selectedFeatures[i]) {
				newAlpha[j] = this.alpha.getEntry(i);
				newBeta[j] = this.beta.getEntry(i);
				newMean[j] = this.mean.getEntry(i);
				newSd[j] = this.sd.getEntry(i);
				newNaRatio[j] = this.naRatio.getEntry(i);
				j++;
			}
		}
		newModel.setAlpha(new ArrayRealVector(Arrays.copyOfRange(newAlpha, 0, j)));
		newModel.setBeta(new ArrayRealVector(Arrays.copyOfRange(newBeta, 0, j)));
		newModel.setMean(new ArrayRealVector(Arrays.copyOfRange(newMean, 0, j)));
		newModel.setSd(new ArrayRealVector(Arrays.copyOfRange(newSd, 0, j)));
		newModel.setNaRatio(new ArrayRealVector(Arrays.copyOfRange(newNaRatio, 0, j)));
		newModel.setnFeature(j);
		newModel.setnSample(nSample);
		return newModel;
	}

	public int getnSample() {
		return nSample;
	}

	public void setnSample(int nSample) {
		this.nSample = nSample;
	}

	public int getnFeature() {
		return nFeature;
	}

	public void setnFeature(int nFeature) {
		this.nFeature = nFeature;
	}

	public RealVector getAlpha() {
		return alpha;
	}

	public RealVector getBeta() {
		return beta;
	}

	public RealVector getMean() {
		return mean;
	}

	public RealVector getSd() {
		return sd;
	}

	public RealVector getNaRatio() {
		return naRatio;
	}


	public void setAlpha(RealVector alpha) { this.alpha = alpha;}

	public void setBeta(RealVector beta) { this.beta = beta;}

	public void setMean(RealVector mean) {
		this.mean = mean;
	}

	public void setSd(RealVector sd) {
		this.sd = sd;
	}

	public void setNaRatio(RealVector naRatio) {
		this.naRatio = naRatio;
	}

	public static double[] estBetaDist(double[] betaValues) {
		Mean mean = new Mean();
		double mu = mean.evaluate(betaValues,0,betaValues.length);
		Variance variance = new Variance();
		double var = variance.evaluate(betaValues, mu);
		double alpha = -mu*(var+mu*mu-mu)/var;
		double beta = (mu-1)*(var+mu*mu-mu)/var;
		return new double[] {alpha, beta, mu, FastMath.sqrt(var)};
	}


}
