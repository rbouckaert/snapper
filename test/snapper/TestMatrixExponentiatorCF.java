package snapper;

import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;

import beast.app.beauti.BeautiDoc;
import junit.framework.TestCase;
import snap.likelihood.COMPLEX;

public class TestMatrixExponentiatorCF extends TestCase {

	COMPLEX [] readComplexArray(String fileName) throws IOException {
		String s = BeautiDoc.load(fileName);
		String [] strs = s.split("\n");
		COMPLEX [] array = new COMPLEX[strs.length];
		for (int i = 0; i < array.length; i++) {
			array[i] = parseComplex(strs[i]);
		}
		return array;
	}


	COMPLEX []readComplexMatrix(String fileName) throws IOException {
		String s = BeautiDoc.load(fileName);
		String [] strs = s.split("\n");
		COMPLEX [][] array = new COMPLEX[strs.length][];
		for (int i = 0; i < array.length; i++) {
			s = strs[i];
			String [] strs2 = s.split("\\)");
			array[i] = new COMPLEX[strs2.length];
			for (int j = 0; j < array[i].length; j++) {
				array[i][j] = parseComplex(strs2[j]);
			}
		}
		
		COMPLEX [] AA = new COMPLEX[array.length * array.length];
		for (int i = 0; i < array.length; i++) {			
			for (int j = 0; j < array.length; j++) {
				AA[i*array.length + j] = array[i][j];
			}
		}
		return AA;
	}
	
	double []readDoubleMatrix(String fileName) throws IOException {
		String s = BeautiDoc.load(fileName);
		String [] strs = s.split("\n");
		int n = strs.length;
		double [] array = new double[n*n];
		for (int i = 0; i < n; i++) {
			s = strs[i];
			String [] strs2 = s.trim().split("\\s+");
			for (int j = 0; j < n; j++) {
				array[i * n + j] = Double.parseDouble(strs2[j]);
			}
		}
		return array;
	}
	
	private COMPLEX parseComplex(String string) {
		int startRe = string.indexOf('(') + 1;
		int endRe = string.indexOf('e') + 3;
		while (string.charAt(endRe) != '+') {
			endRe++;
		}
		int startIm = endRe + 1;
		int endIm = string.indexOf('j');
		String re = string.substring(startRe,endRe);
		String im = string.substring(startIm,endIm);
		try {
			return new COMPLEX(Double.parseDouble(re), Double.parseDouble(im));
		} catch (NumberFormatException e) {
			return null;
		}
	}

	@Test
	public void testGetValuesTridag() throws IOException {
		COMPLEX [] AA = readComplexMatrix("test/snapper/test1/AA.txt");

		COMPLEX [] A2c = readComplexArray("test/snapper/test1/A2c.txt");
		COMPLEX [] a = readComplexArray("test/snapper/test1/a.txt");
		COMPLEX [] b = readComplexArray("test/snapper/test1/b.txt");
		COMPLEX [] c = readComplexArray("test/snapper/test1/c.txt");
		COMPLEX [] d = readComplexArray("test/snapper/test1/d.txt");
		
		MatrixExponentiator e = new MatrixExponentiator();
		COMPLEX [][] out = new COMPLEX[4][];
		out[0] = new COMPLEX[a.length]; init(out[0]);
		out[1] = new COMPLEX[b.length]; init(out[1]);
		out[2] = new COMPLEX[c.length]; init(out[2]);
		out[3] = new COMPLEX[d.length]; init(out[3]);
		
		e.getValuesTridag(AA, A2c, out);
		
		compare(out[0], a);
		compare(out[1], b);
		compare(out[2], c);
		compare(out[3], d);
	}

	@Test
	public void testGetValuesUBand() throws IOException {
		COMPLEX [] AA = readComplexMatrix("test/snapper/test1/AA.txt");
		double [] LL = readDoubleMatrix("test/snapper/test1/LL.txt");

		COMPLEX [] c_i = readComplexArray("test/snapper/test1/c_i.txt");

		COMPLEX [] A2c = readComplexArray("test/snapper/test1/A2c.txt");
		COMPLEX [] a = readComplexArray("test/snapper/test1/a0.txt");
		COMPLEX [] b = readComplexArray("test/snapper/test1/b0.txt");
		COMPLEX [] c = readComplexArray("test/snapper/test1/c0.txt");
		COMPLEX [] d = readComplexArray("test/snapper/test1/d0.txt");
		
		MatrixExponentiator e = new MatrixExponentiator();
		COMPLEX [][] out = new COMPLEX[4][];
		out[0] = new COMPLEX[a.length]; init(out[0]);
		out[1] = new COMPLEX[b.length]; init(out[1]);
		out[2] = new COMPLEX[c.length]; init(out[2]);
		out[3] = new COMPLEX[d.length]; init(out[3]);
		
		COMPLEX z = new COMPLEX(-699.8688082445777,1399.5917029301354);
	
		
		e.getValuesUBand(LL, AA, A2c, c_i, z, out);
		
		compare(out[0], a);
		compare(out[1], b);
		compare(out[2], c);
		compare(out[3], d);
	}

	@Test
	public void testUBandSolver() throws IOException {
		COMPLEX [] a = readComplexArray("test/snapper/test1/a0.txt");
		COMPLEX [] b = readComplexArray("test/snapper/test1/b0.txt");
		COMPLEX [] c = readComplexArray("test/snapper/test1/c0.txt");
		COMPLEX [] d = readComplexArray("test/snapper/test1/d0.txt");
		
		COMPLEX [] solvEven = readComplexArray("test/snapper/test1/sol_even.txt");
		COMPLEX [] x = new COMPLEX[solvEven.length]; init(x);
		
		MatrixExponentiator e = new MatrixExponentiator();
		COMPLEX [][] in = new COMPLEX[4][];
		in[0] = a;
		in[1] = b;
		in[2] = c;
		in[3] = d;

		e.ubandSolver(in, x, 33);
		
		compare(x, solvEven);
	}
	
	@Test
	public void testTriDagSolver() throws IOException {
		COMPLEX [] a = readComplexArray("test/snapper/test1/a.txt");
		COMPLEX [] b = readComplexArray("test/snapper/test1/b.txt");
		COMPLEX [] c = readComplexArray("test/snapper/test1/c.txt");
		COMPLEX [] d = readComplexArray("test/snapper/test1/d.txt");
		
		COMPLEX [] solvOdd = readComplexArray("test/snapper/test1/sol_odd.txt");
		COMPLEX [] x = new COMPLEX[solvOdd.length]; init(x);
		
		MatrixExponentiator e = new MatrixExponentiator();
		COMPLEX [][] in = new COMPLEX[4][];
		in[0] = a;
		in[1] = b;
		in[2] = c;
		in[3] = d;

		e.triDiagSolver(in, x, 33);
		
		compare(x, solvOdd);
	}

	@Test
	public void testFasterSolver() throws IOException {
		COMPLEX [] AA = readComplexMatrix("test/snapper/test1/AA.txt");
		double [] LL = readDoubleMatrix("test/snapper/test1/LL.txt");

		COMPLEX [] c_i = readComplexArray("test/snapper/test1/c_i.txt");

		COMPLEX [] A2c = readComplexArray("test/snapper/test1/A2c.txt");
		
		COMPLEX z = new COMPLEX(-699.8688082445777,1399.5917029301354);

		MatrixExponentiator e = new MatrixExponentiator();

		COMPLEX [] x = e.fasterSolver(LL, AA, A2c, c_i, z);
		
		COMPLEX [] solv = readComplexArray("test/snapper/test1/solve.txt");
		compare(x, solv);
	}

	@Test
	public void testExpCF() throws IOException {

		MatrixExponentiator e = new MatrixExponentiator();
		double [] v = new double[33];
		QMatrix Q = new QMatrix(33);
		double tt = 1;
		double u_ = 0.001;
		double v_ = 0.001;
		double [] a = new double[]{u_,-(u_ + v_),0,0};
		double [] b = new double[]{0,1,-1,0};

		Q.setQ(a, b);
		e.expCF(v, Q.Q, tt);
		double [] solv = new double[] {1.52243470e+001, 1.58564913e+001, 1.65842228e+000, 2.44590626e-002, 5.44686091e-005, 1.58507873e-008, 6.76217216e-013, 3.70858069e-018, 2.95950235e-024, 3.01419867e-031, 4.45878137e-039, 8.37584314e-048, 2.28665015e-057, 7.89072733e-068, 3.96602854e-079, 2.50676180e-091, 2.31619252e-104, 2.67464499e-118, 4.53870264e-133, 9.54988368e-149, 2.97428242e-165, 1.13656637e-182, 6.49374835e-201, 4.48569695e-220, 4.70002123e-240, 5.82483864e-261, 1.11895592e-282, 2.45305418e-305, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000};
		compare(v, solv);
	}
	

	
	private void compare(COMPLEX[] out, COMPLEX[] d) {
		assertEquals(out.length, d.length);
		for (int i = 0; i < out.length; i++) {
			assertEquals(d[i].m_fRe, out[i].m_fRe, 1e-10);
			assertEquals(d[i].m_fIm, out[i].m_fIm, 1e-10);
		}
	}

	private void compare(double[] out, double[] d) {
		assertEquals(out.length, d.length);
		for (int i = 0; i < out.length; i++) {
			assertEquals(d[i], out[i], 1e-10);
		}
	}

	private void init(COMPLEX[] array) {
		for (int i = 0; i < array.length; i++) {
			array[i] = new COMPLEX();
		}
	}
	
	
	
}
