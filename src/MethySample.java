import java.util.Arrays;

public class MethySample {
	private String id;
	private int index;
	private String type;
	private boolean realData;
	private double theta;
	private int[] depth;
	private double[] methy;
	// for predicted samples
	private  double densRatio; // likelihood ratio

	public MethySample(String id, String type, int index, double theta, boolean realData) {
		this.id = id;
	    this.type = type;
	    this.index = index;
	    this.theta = theta;
	    this.realData = realData;
	}
	
	public MethySample(String type, double theta) {
	    this.type = type;
	    this.theta = theta;
	}
	
	public MethySample() {
	}

	public void selfSelectFeature(boolean[] selectedFeatures) {
		int maxFeature = selectedFeatures.length;
		int[] newDepth = new int[maxFeature];
		double[] newMethy = new double[maxFeature];
		int j = 0; // num of selected features
		for(int i=0; i<maxFeature; i++) {
			if (selectedFeatures[i]) {
				newDepth[j] = this.depth[i];
				newMethy[j] = this.methy[i];
				j++;
			}
		}
		this.setDepth(Arrays.copyOfRange(newDepth, 0, j));
		this.setMethy(Arrays.copyOfRange(newMethy, 0, j));
	}

	public MethySample selectFeature(boolean[] selectedFeatures) {
		MethySample newSample = new MethySample(this.id, this.type,this.index,this.theta, this.realData);
		int maxFeature = selectedFeatures.length;
		int[] newDepth = new int[maxFeature];
		double[] newMethy = new double[maxFeature];
		int j = 0; // num of selected features
		for(int i=0; i<maxFeature; i++) {
			if (selectedFeatures[i]) {
				newMethy[j] = this.methy[i];
				newDepth[j] = this.depth[i];
				j++;
			}
		}
		newSample.setDepth(Arrays.copyOfRange(newDepth, 0, j));
		newSample.setMethy(Arrays.copyOfRange(newMethy, 0, j));
		return newSample;
	}

	public MethySample selectFeature(int[] selectedFeatures) {
		MethySample newSample = new MethySample(this.id, this.type,this.index,this.theta, this.realData);
		int nFeature = selectedFeatures.length;
		double[] newMethy = new double[nFeature];
		int[] newDepth = new int[nFeature];
		for(int i=0; i<nFeature; i++) {
			newMethy[i] = this.methy[selectedFeatures[i]];
			newDepth[i] = this.depth[selectedFeatures[i]];
		}
		newSample.setMethy(newMethy);
		newSample.setDepth(newDepth);
		return newSample;
	}


	public int getFeatureNum() {
		return this.methy.length;
	}
	
	public void setType(String type) {
		this.type = type;
	}
	public String getType() {
		return type;
	}
	
	public void setIndex(int index) {
		this.index = index;
	}
	public int getIndex() {
		return index;
	}

	public void setTheta(double theta) {
		this.theta = theta;
	}
	public double getTheta() {
		return theta;
	}

	public boolean isRealData() {
		return realData;
	}

	public void setDepth(int[] depth) {this.depth = depth;}
	public int[] getDepth() {return depth;}

	public void setMethy(double[] methy) {
		this.methy = methy;
	}
	public double[] getMethy() {
		return methy;
	}

	public String getId() {
		return id;
	}

	public double getDensRatio() {
		return densRatio;
	}

	public void setDensRatio(double densRatio) {
		this.densRatio = densRatio;
	}

	public void setId(String id) {
		this.id = id;
	}

	public void setRealData(boolean realData) {
		this.realData = realData;
	}
}
