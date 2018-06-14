package snapper;

import java.util.Arrays;

import org.apache.commons.math3.transform.DctNormalization;
import org.apache.commons.math3.transform.FastCosineTransformer;
import org.apache.commons.math3.transform.TransformType;

public class ChebyshevPolynomial {	
	static FastCosineTransformer transformer = null;
	
	double [] f; // function values
	double [] a; // Chebyshev polynomial coefficients

	public ChebyshevPolynomial(int N) {
		f = new double[N];
		a = new double[N];
		if (transformer == null) {
			transformer = new org.apache.commons.math3.transform.FastCosineTransformer(DctNormalization.STANDARD_DCT_I);
		}
	}
	
	/**
	 * initialises f and a to approximate a binomial(r,n) distribution 
	 * @param r number of red lineages
	 * @param n total number of lineages
	 */
	public void init(int r, int n) {
		double c = logBinom(r, n);
		int N = f.length;
//		double [] f = new double[n]; 
		int nf = f.length;


		
		if (n == 0) { 
			// deal with missing data
			Arrays.fill(a, 0);
			a[0] = 1;
			aToF();
		} else {
			// there is some data
	    	// first, set f[0] and f[N-1];
	    	if (r == 0) {
	    		f[0] = 1;
	    		f[nf-1] = 0;            		
	    	} else if (r == n) {
	    		f[0] = 0;
	    		f[nf-1] = 1;
	    	} else {
	    		f[0] = 0;
	    		f[nf-1] = 0;
	    	}
	    		
	    	// set the other values
	    	for (int m = 1; m < nf - 1; m++) {
	    		double x = 0.5 - Math.cos(-m/(nf-1.0)*Math.PI) / 2.0;
	    		double logp = c + r * Math.log(x) + (n-r) * Math.log(1-x);
	    		double p = Math.exp(logp);
	    		f[m] = p;
	    	}
	    	fToA();
//	    	// slow discrete cosine transform
//	    	for (int i = 0; i < N; i++) {
//	    		double sum = 0;
//	    		for (int j = 0; j < n; j++) {
//	    			sum += f[j] * Math.cos((2*j+1)*i*Math.PI/(2*N));
//	    		}
//	    		if (i == 0) {
//	    			sum *= 1.0/n;//Math.sqrt(n);
//	    		} else {
//	    			sum *= 2.0/n;//Math.sqrt(n);
//	    		}
//	    		a[i] = sum;
//	    	}
//	    	System.out.println(Arrays.toString(f));
//			aToF();
		}
	}
	
	// convert coefficients to function values 
	public void aToF() {
		a[0] *= 2.0;
		a[a.length - 1] *= 2.0;
		f = transformer.transform(a, TransformType.FORWARD);
		a[0] /= 2.0;
		a[a.length - 1] /= 2.0;
	}

	// convert function values to coefficients
	public void fToA() {		
		a = transformer.transform(f, TransformType.INVERSE);
		a[0] /= 2.0;
		a[a.length - 1] /= 2.0;
 	}

	public void setPolyFactors(double[] a) {
		if (this.a.length != a.length) {
			throw new IllegalArgumentException("Expected dimension " + this.a.length + " but got " + a.length);
		}
		System.arraycopy(a, 0, this.a, 0, a.length);
	}

	public void setPolyValues(double[] f) {
		if (this.f.length != f.length) {
			throw new IllegalArgumentException("Expected dimension " + this.f.length + " but got " + f.length);
		}
		System.arraycopy(f, 0, this.f, 0, f.length);
	}
	
    private double logBinom(int k, int n) {
    	double f = 0;
    	for (int i = k + 1; i <= n; i++) {
    		f += Math.log(i) - Math.log(n - i + 1);
    	}
		return f;
	}

    @Override
    public String toString() {
    	return "a: " + Arrays.toString(a) + "\n" +
    		   "f: " + Arrays.toString(f) + "\n";
    }
}
