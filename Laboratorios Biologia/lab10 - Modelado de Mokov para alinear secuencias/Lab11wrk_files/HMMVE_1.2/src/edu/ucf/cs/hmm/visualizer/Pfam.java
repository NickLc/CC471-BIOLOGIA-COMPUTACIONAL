package edu.ucf.cs.hmm.visualizer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;

class PfamIndexEntry {
	String hmmName;
	String hmmDesc;
	long startFilePointer;
}

interface IndexingListener {
	void indexProgress(int i);
	void indexFinish();
}

/**
 * 
 * Structure to store Pfam index
 *
 */
public class Pfam {
	String pfamFileName;
	boolean complete = false;
	long fileSize;
	long position;
	boolean canceled;
	public Vector<PfamIndexEntry> indexEntries = new Vector<PfamIndexEntry>();
	IndexingListener listener = null;
	
	public Pfam(String f) throws IOException
	{
		File pfamFile = new File(f);
		if (!pfamFile.canRead())
			throw new IOException();
		pfamFileName = f;
		fileSize = pfamFile.length();
	}
	public void addIndexingListener(IndexingListener l)
	{
		listener = l;
	}
	public void buildIndex() throws IOException
	{
		complete = false;
		canceled = false;
		BufferedReader reader =  new BufferedReader(new FileReader(pfamFileName));
		String buffer;
		position = 0;
		PfamIndexEntry newEntry=null;
		while ((buffer=reader.readLine())!=null && !canceled)
		{
			if (buffer.startsWith("HMMER"))
			{
				newEntry = new PfamIndexEntry();
				newEntry.startFilePointer=position;
			}
			if (buffer.startsWith("NAME"))  // name of current entry
			{
				if (newEntry.hmmName==null)  // Always
				{
					newEntry.hmmName=buffer.split("\\s+")[1];
				}
			}
			if (buffer.startsWith("DESC"))
			{
				if (newEntry.hmmDesc==null)  // Always
				{
					newEntry.hmmDesc=buffer.substring(6);
				}
			}
			if (buffer.startsWith("//"))  // End of hmm
			{
				indexEntries.add(newEntry);
				if (indexEntries.size()%100==0 && listener!=null)
					listener.indexProgress(indexEntries.size());
			}
			position = position + buffer.length()+1;
		}
		reader.close();
		if (!canceled)
		{
			if (listener!=null)
				listener.indexFinish();
			complete = true;
		}
	}
	void reset()
	{
		indexEntries.clear();
		complete = false;
		position = 0;
	}
	void cancel()
	{
		canceled = true;
	}
	Hmm getHmm(int index) throws IOException, InvalidFormatException
	{
		Hmm result = new Hmm();
		BufferedReader reader = new BufferedReader(new FileReader(pfamFileName));
		reader.skip(indexEntries.get(index).startFilePointer);
		result.read(reader);
		reader.close();
		return result;
	}
}

