package edu.ucf.cs.hmm.squid;

import java.io.BufferedReader;

public class SSIFILE {
	BufferedReader reader; /* open SSI index file                 */

	int flags; /* optional behavior flags             */

	int nfiles; /* number of files = 16 bit int        */

	int nprimary; /* number of primary keys              */

	int nsecondary; /* number of secondary keys            */

	int flen; /* length of filenames (inc '\0')      */

	int plen; /* length of primary keys (inc '\0')   */

	int slen; /* length of secondary keys (inc '\0') */

	int frecsize; /* # bytes in a file record            */

	int precsize; /* # bytes in a primary key record     */

	int srecsize; /* # bytes in a secondary key record   */

	int foffset; /* disk offset, start of file records  */

	int poffset; /* disk offset, start of pri key recs  */

	int soffset; /* disk offset, start of sec key recs  */

	char imode; /* mode for index file offsets, 32 v. 64 bit    */

	char smode; /* mode for sequence file offsets, 32 v. 64 bit */

	/* File information:
	 */
	String[] filename; /* list of file names [0..nfiles-1]    */

	int fileformat; /* file formats                        */

	int fileflags; /* optional per-file behavior flags    */

	int bpl; /* bytes per line in file              */

	int rpl; /* residues per line in file           */
}
