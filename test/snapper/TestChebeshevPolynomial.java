package snapper;

import org.junit.Test;

import junit.framework.TestCase;

public class TestChebeshevPolynomial extends TestCase {
	
	@Test
	public void testChebeshevPolynomial() {
		int N = 33;
		ChebyshevPolynomial p = new ChebyshevPolynomial(N);
		p.init(60, 100);
		
		// store function values
		double [] f = new double[N];
		System.arraycopy(p.f, 0, f, 0, f.length);
		
		// convert factors to function values
		p.aToF();
		
		// check for largest error in function
		double maxErr = 0;
		int iErr = -1;
		for (int i = 0; i < f.length; i++) {
			if (Math.abs(f[i] - p.f[i]) > maxErr) {
				maxErr = Math.abs(f[i] - p.f[i]);
				iErr = i;
			}
		}
		System.out.println("Max error: " + maxErr + "\nat entry " + iErr);
		assertTrue(Math.abs(maxErr) < 1e-14);
	}

}
