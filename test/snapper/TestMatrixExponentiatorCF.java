package snapper;

import java.io.IOException;

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

	private void compare(COMPLEX[] out, COMPLEX[] d) {
		assertEquals(out.length, d.length);
		for (int i = 0; i < out.length; i++) {
			assertEquals(d[i].m_fRe, out[i].m_fRe, 1e-10);
			assertEquals(d[i].m_fIm, out[i].m_fIm, 1e-10);
		}
	}


	private void init(COMPLEX[] array) {
		for (int i = 0; i < array.length; i++) {
			array[i] = new COMPLEX();
		}
	}
	
	
	
}
