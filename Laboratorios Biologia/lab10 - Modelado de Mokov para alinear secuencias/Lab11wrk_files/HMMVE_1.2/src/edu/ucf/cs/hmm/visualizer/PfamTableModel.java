package edu.ucf.cs.hmm.visualizer;

import java.util.Vector;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import javax.swing.event.TableModelListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

/**
 * Table Model for pfam index
 * 
 */
class PfamTableModel extends AbstractTableModel implements TableModel
{
	Pfam pfam;
	String filter;
	int filterSoFar;
	Vector<TableModelListener> listeners = new Vector<TableModelListener>();
	Vector<Integer> filtered = new Vector<Integer>();
	public PfamTableModel(Pfam p)
	{
		pfam = p;
		filter = null;
		filterSoFar = -1;
	}
	public void setFilter(String f)
	{
		if ("".equals(f))
		{
			filter = null;
			filtered.clear();
			return;
		}
		Pattern filterPattern=null;
		try {
			filterPattern = Pattern.compile(f, Pattern.CASE_INSENSITIVE);
			filter = f;
		} catch (PatternSyntaxException e) {
			return;
		}
		synchronized(filtered)
		{
			filtered.clear();
			filterSoFar = pfam.indexEntries.size();
			for (int i=0;i<filterSoFar;i++)
			{
				if (filterPattern.matcher(pfam.indexEntries.get(i).hmmName).find() ||
						filterPattern.matcher(pfam.indexEntries.get(i).hmmDesc).find())
				{
					filtered.add(i);
				}
			}
		}
	}
	public void processNewEntries()
	{
		if (filter!=null && filterSoFar!=-1)  // let batch process first
		{
			Pattern filterPattern = Pattern.compile(filter);
			synchronized(filtered)
			{
				int low = filterSoFar;
				filterSoFar = pfam.indexEntries.size();
				for (int i=low;i<filterSoFar;i++)
				{
					if (filterPattern.matcher(pfam.indexEntries.get(i).hmmName).find() ||
							filterPattern.matcher(pfam.indexEntries.get(i).hmmDesc).find())
					{
						filtered.add(i);
					}
				}
			}	
		}
	}
	public void reset()
	{
		filtered.clear();
		filter = null;
		filterSoFar = -1;
	}
	public void addTableModelListener(TableModelListener l) 
	{
		listeners.add(l);
	}
	public Class getColumnClass(int columnIndex)
	{
		return String.class;
	}
	public int getColumnCount()
	{
		return 3;
	}
	public String getColumnName(int columnIndex) 
	{
		if (columnIndex==0)
			return "#";
		if (columnIndex==1)
			return "Hmm Name";
		if (columnIndex==2)
			return "Description";
		return "Opps, Wrong entry!";
	}
	public int getRowCount()
	{
		if (filter==null)
			return pfam.indexEntries.size();
		return filtered.size();
	}
	public Object getValueAt(int rowIndex, int columnIndex)
	{
		if (filter==null)
		{
			if (columnIndex==0)
				return rowIndex;
			else if (columnIndex==1)
				return pfam.indexEntries.get(rowIndex).hmmName;
			else if (columnIndex==2)
				return pfam.indexEntries.get(rowIndex).hmmDesc;
			else
				return "Opps, what is it?";
		}
		//filter is not null
		if (columnIndex==0)
			return filtered.get(rowIndex);
		else if (columnIndex==1)
			return pfam.indexEntries.get(filtered.get(rowIndex)).hmmName;
		else if (columnIndex==2)
			return pfam.indexEntries.get(filtered.get(rowIndex)).hmmDesc;
		else
			return "Opps, what happens?";
	}
	public boolean isCellEditable(int rowIndex, int columnIndex) 
	{
		return false;
	}
	public void removeTableModelListener(TableModelListener l)
	{
		listeners.remove(l);
	}
	public void setValueAt(Object aValue, int rowIndex, int columnIndex)
	{
	}
}
