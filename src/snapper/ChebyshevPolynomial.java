package snapper;

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
			transformer = new FastCosineTransformer(DctNormalization.STANDARD_DCT_I);
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
}
