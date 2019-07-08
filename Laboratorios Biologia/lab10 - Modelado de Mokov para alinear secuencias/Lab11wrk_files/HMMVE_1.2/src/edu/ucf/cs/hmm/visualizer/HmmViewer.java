package edu.ucf.cs.hmm.visualizer;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import javax.imageio.ImageIO;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.ScrollPaneConstants;
import javax.swing.border.EtchedBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;

import edu.ucf.cs.hmm.squid.MSA;
import edu.ucf.cs.hmm.squid.MSAFILE;
import edu.ucf.cs.hmm.squid.SQFILE;
import edu.ucf.cs.hmm.squid.Sequence;

/**
 * Main frame of HMMVE
 */
public class HmmViewer extends JFrame implements ActionListener,
		ComponentListener, ChangeListener, WindowListener, MouseListener {

	private static final long serialVersionUID = 1L;

	private JMenuBar jJMenuBar = null;

	private JSplitPane splitPane = null;

	private JPanel lowerPane = null;

	private JScrollPane jHmmScrollPane = null;

	private JScrollPane jNullScrollPane = null;

	private JScrollPane jHmmTextScrollPane = null;

	private JScrollPane jLogoScrollPane = null;

	JScrollPane sequencesScrollpane = null;

	JTable alignResultTable = null;

	JScrollPane alignResultScrollPane = null;

	private HmmPanel hmmCanvas = null;

	private NullModelPanel nullCanvas = null;

	private LogoPanel logoCanvas = null;

	JScrollPane hmmTextScrollPane = null;

	HmmText hmmText = new HmmText();

	private JTabbedPane tabbedHmmPane = null;

	private JTabbedPane tabbedInfoPane = null;

	private JPanel sequencePane = null;

	private JPanel resultPane = null;

	private StateChartLabel jChartLabel = null;

	private JTable jSequenceTable = null;

	private JButton jAlignPathButton = null;

	private JButton jAlignAllButton = null;

	private JButton jClearButton = null;

	private JMenuItem undoMenuItem = null; // need this handle to enable and

	// disable it

	private JMenuItem redoMenuItem = null; // need this handle to enable and

	// disable it

	private JCheckBoxMenuItem hidePathMenuItem = null; // will check the

	// checking state of
	// this menu item

	private JCheckBoxMenuItem hideEmissionMenuItem = null; // will check the

	// checking state of
	// this menu item

	private JButton jSaveSequenceButton = null;

	private String hmmFileName;
	
	private String pfamFileName;

	private Hmm hmm = null;
	
	private Pfam pfam = null;

	private String hmmString;

	private HmmLayout hmmLayout;

	private NullModelLayout nullLayout;

	private HmmCommandList commandHistory;

	MSA msa = null;

	static final double ZoomFactor = 1.2;

	int precision = 2;

	static final int MinPrecision = 1;

	static final int MaxPrecision = 10;

	private String sequenceFileName;

	private ArrayList<Sequence> sequences = null;

    private File currentDirectory=  new File(".");

	class HmmFileFilter extends FileFilter {
		public boolean accept(File f) {
			if (f.isDirectory()) {
				return true;
			} else if (f.getName().toUpperCase().endsWith("HMM"))
				return true;
			return false;
		}

		public String getDescription() {
			return "Hmm File";
		}
	};

	class SequenceFileFilter extends FileFilter {
		public boolean accept(File f) {
			if (f.isDirectory()) {
				return true;
			} else if (f.getName().toUpperCase().endsWith("MSF")
					|| f.getName().toUpperCase().endsWith("FA")
					|| f.getName().toUpperCase().endsWith("STO")
					|| f.getName().toUpperCase().endsWith("TFA"))
				return true;
			return false;
		}

		public String getDescription() {
			return "All Supported Sequence File(*.msf, *.fa, *.sto, *.tfa)";
		}
	};

	class MSFFileFilter extends FileFilter {
		public boolean accept(File f) {
			if (f.isDirectory()) {
				return true;
			} else if (f.getName().toUpperCase().endsWith("MSF"))
				return true;
			return false;
		}

		public String getDescription() {
			return "MSF Sequence File(*.msf)";
		}
	}

	/**
	 * This is the default constructor
	 */
	public HmmViewer() {
		super();
		initialize();
	}

	/**
	 * This method initializes this
	 * 
	 * @return void
	 */
	private void initialize() {
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		this.setSize(800, 600);
		this.setResizable(true);
		this.setTitle("HMMVE");
		this.setJMenuBar(getJJMenuBar());
		commandHistory = new HmmCommandList(undoMenuItem, redoMenuItem);
		this.setContentPane(getJContentPane());
		this.addComponentListener(this);
		resize();
		this.addWindowListener(this);
	}

	/**
	 * This method initializes jContentPane
	 * 
	 * @return javax.swing.JPanel
	 */
	private JSplitPane getJContentPane() {
		if (splitPane == null) {
			jHmmScrollPane = new JScrollPane(
					ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
					ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
			jHmmScrollPane.getViewport().setBackground(Color.WHITE);
			jHmmScrollPane.setWheelScrollingEnabled(true);

			jNullScrollPane = new JScrollPane(
					ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
					ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
			jNullScrollPane.getViewport().setBackground(Color.WHITE);
			jNullScrollPane.setWheelScrollingEnabled(true);

			jHmmTextScrollPane = new JScrollPane(
					ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
					ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
			jHmmTextScrollPane.setWheelScrollingEnabled(true);
			hmmTextScrollPane = new JScrollPane(hmmText,
					ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
					ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);

			jLogoScrollPane = new JScrollPane(
					ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
					ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
			jLogoScrollPane.getViewport().setBackground(Color.WHITE);
			jLogoScrollPane.setWheelScrollingEnabled(true);

			tabbedHmmPane = new JTabbedPane();
			tabbedHmmPane.addTab("Hmm Model", jHmmScrollPane);
			tabbedHmmPane.addTab("Null Model", jNullScrollPane);
			tabbedHmmPane.addTab("Hmm Logo", jLogoScrollPane);
			tabbedHmmPane.addTab("Hmm Text", hmmTextScrollPane);
			tabbedHmmPane.addChangeListener(this);

			jChartLabel = new StateChartLabel();
			jChartLabel.setOpaque(true);
			jChartLabel.setBorder(BorderFactory
					.createEtchedBorder(EtchedBorder.RAISED));
			jChartLabel.setBackground(Color.WHITE);

			jSequenceTable = new JTable();
			jSequenceTable.addMouseListener(this);
			sequencesScrollpane = new JScrollPane(jSequenceTable,
					JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
					JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

			jAlignPathButton = new JButton("Align Path");
			jAlignPathButton.setActionCommand("AlignPath");
			jAlignPathButton.addActionListener(this);
			jAlignPathButton.setBorder(BorderFactory.createRaisedBevelBorder());

			jClearButton = new JButton("Clear Path");
			jClearButton.setActionCommand("Clear");
			jClearButton.addActionListener(this);
			jClearButton.setBorder(BorderFactory.createRaisedBevelBorder());

			jAlignAllButton = new JButton("Align All");
			jAlignAllButton.setActionCommand("AlignAll");
			jAlignAllButton.addActionListener(this);
			jAlignAllButton.setBorder(BorderFactory.createRaisedBevelBorder());

			sequencePane = new JPanel();
			sequencePane.setLayout(null);

			sequencePane.add(sequencesScrollpane);
			sequencePane.add(jAlignPathButton);
			sequencePane.add(jClearButton);
			sequencePane.add(jAlignAllButton);

			alignResultTable = new JTable();
			alignResultTable.setFont(new Font("Courier New", Font.BOLD, 12));
			alignResultTable.addMouseListener(this);
			alignResultScrollPane = new JScrollPane(alignResultTable,
					JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
					JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

			jSaveSequenceButton = new JButton("Save Sequence");
			jSaveSequenceButton.setActionCommand("SaveSequence");
			jSaveSequenceButton.addActionListener(this);
			jSaveSequenceButton.setBorder(BorderFactory
					.createRaisedBevelBorder());

			resultPane = new JPanel();
			resultPane.setLayout(null);
			resultPane.add(alignResultScrollPane);
			resultPane.add(jSaveSequenceButton);

			tabbedInfoPane = new JTabbedPane();
			tabbedInfoPane.addTab("Sequences", sequencePane);
			tabbedInfoPane.addTab("Result", resultPane);

			lowerPane = new JPanel();
			lowerPane.setLayout(null);
			lowerPane.add(jChartLabel);
			lowerPane.add(tabbedInfoPane);
			lowerPane.addComponentListener(this);

			splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
					tabbedHmmPane, lowerPane);
			splitPane.setResizeWeight(0.6);
		}
		return splitPane;
	}

	/**
	 * This method initializes jJMenuBar
	 * 
	 * @return javax.swing.JMenuBar
	 */
	private JMenuBar getJJMenuBar() {
		if (jJMenuBar == null) {
			jJMenuBar = new JMenuBar();
		}
		JMenu menu = new JMenu("File");
		menu.setMnemonic(KeyEvent.VK_F);
		JMenuItem menuItem = new JMenuItem("Open Hmm...", KeyEvent.VK_O);
		menuItem.setAccelerator(KeyStroke.getKeyStroke("ctrl O"));
		menuItem.setActionCommand("OpenHmm");
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menuItem = new JMenuItem("Open Pfam...", KeyEvent.VK_P);
		menuItem.setAccelerator(KeyStroke.getKeyStroke("ctrl P"));
		menuItem.setActionCommand("OpenPfam");
		menuItem.addActionListener(this);
		menu.add(menuItem);
		
		menuItem = new JMenuItem("Save Hmm...", KeyEvent.VK_S);
		menuItem.setActionCommand("SaveHmm");
		menuItem.setAccelerator(KeyStroke.getKeyStroke("ctrl S"));
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menuItem = new JMenuItem("Open Sequence...", KeyEvent.VK_E);
		menuItem.setActionCommand("OpenSequence");
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menuItem = new JMenuItem("Save Image...", KeyEvent.VK_I);
		menuItem.setActionCommand("SaveImage");
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menuItem = new JMenuItem("Exit", KeyEvent.VK_E);
		menuItem.setActionCommand("Exit");
		menuItem.addActionListener(this);
		menu.add(menuItem);
		jJMenuBar.add(menu);

		menu = new JMenu("Edit");
		menu.setMnemonic(KeyEvent.VK_E);
		undoMenuItem = new JMenuItem("Undo", KeyEvent.VK_U);
		undoMenuItem.setAccelerator(KeyStroke.getKeyStroke("ctrl Z"));
		undoMenuItem.setActionCommand("Undo");
		undoMenuItem.addActionListener(this);
		undoMenuItem.setEnabled(false);
		menu.add(undoMenuItem);

		redoMenuItem = new JMenuItem("Redo", KeyEvent.VK_R);
		redoMenuItem.setAccelerator(KeyStroke.getKeyStroke("ctrl Y"));
		redoMenuItem.setActionCommand("Redo");
		redoMenuItem.setEnabled(false);
		redoMenuItem.addActionListener(this);
		menu.add(redoMenuItem);

		jJMenuBar.add(menu);

		menu = new JMenu("View");
		menu.setMnemonic(KeyEvent.VK_V);
		menuItem = new JMenuItem("Zoom in", KeyEvent.VK_I);
		menuItem.setActionCommand("ZoomIn");
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menuItem = new JMenuItem("Zoom out", KeyEvent.VK_O);
		menuItem.setActionCommand("ZoomOut");
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menu.addSeparator();

		menuItem = new JMenuItem("Increase precision", KeyEvent.VK_P);
		menuItem.setActionCommand("IncreasePrecision");
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menuItem = new JMenuItem("Decrease precision", KeyEvent.VK_E);
		menuItem.setActionCommand("DecreasePrecision");
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menu.addSeparator();

		menuItem = new JMenuItem("Default layout", KeyEvent.VK_D);
		menuItem.setActionCommand("DefaultLayout");
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menuItem = new JMenuItem("Re-Number", KeyEvent.VK_R);
		menuItem.setActionCommand("ReNumber");
		menuItem.addActionListener(this);
		menu.add(menuItem);

		menu.addSeparator();

		ButtonGroup group = new ButtonGroup();

		menuItem = new JRadioButtonMenuItem("View label");
		menuItem.setMnemonic(KeyEvent.VK_L);
		menuItem.setActionCommand("ViewLabel");
		menuItem.setSelected(true);
		menuItem.addActionListener(this);
		group.add(menuItem);
		menu.add(menuItem);

		menuItem = new JRadioButtonMenuItem("View consensus");
		menuItem.setMnemonic(KeyEvent.VK_C);
		menuItem.setActionCommand("ViewConsensus");
		menuItem.setSelected(false);
		menuItem.addActionListener(this);
		group.add(menuItem);
		menu.add(menuItem);

		menu.addSeparator();

		hidePathMenuItem = new JCheckBoxMenuItem("Hide redundant path", false);
		hidePathMenuItem.setMnemonic(KeyEvent.VK_H);
		hidePathMenuItem.setActionCommand("HidePath");
		hidePathMenuItem.addActionListener(this);
		menu.add(hidePathMenuItem);

		hideEmissionMenuItem = new JCheckBoxMenuItem("Hide emission string",
				false);
		hideEmissionMenuItem.setMnemonic(KeyEvent.VK_S);
		hideEmissionMenuItem.setActionCommand("HideEmission");
		hideEmissionMenuItem.addActionListener(this);
		menu.add(hideEmissionMenuItem);

		jJMenuBar.add(menu);

		return jJMenuBar;
	}

	/**
	 * When resize the main window
	 *
	 */
	void resize() {
		Rectangle r = lowerPane.getBounds();
		jChartLabel.setBounds(new Rectangle(10, 10, r.height - 20,
				r.height - 20));
		tabbedInfoPane.setBounds(new Rectangle(r.height, 10, r.width - r.height
				- 5, r.height - 20));

		Rectangle s = tabbedInfoPane.getBounds();
		sequencesScrollpane.setBounds(new Rectangle(10, 10, s.width - 110,
				s.height - 40));
		jAlignPathButton.setBounds(new Rectangle(s.width - 95, 10, 90, 20));
		jAlignPathButton.setFont(new Font("Serif", Font.BOLD, 12));
		jClearButton.setBounds(new Rectangle(s.width - 95, 40, 90, 20));
		jClearButton.setFont(new Font("Serif", Font.BOLD, 12));
		jAlignAllButton.setBounds(new Rectangle(s.width - 95, 70, 90, 20));
		jAlignAllButton.setFont(new Font("Serif", Font.BOLD, 12));

		alignResultScrollPane.setBounds(new Rectangle(10, 10, s.width - 110,
				s.height - 40));
		jSaveSequenceButton.setBounds(new Rectangle(s.width - 95, 10, 95, 20));
	}

	public static void main(String[] args) {
		// Schedule a job for the event-dispatching thread:
		// creating and showing this application's GUI.
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				JFrame.setDefaultLookAndFeelDecorated(true);
				HmmViewer viewer = new HmmViewer();
				viewer.setVisible(true);
			}
		});
	}

	public void initWindowLoadHmm()
	{
		tabbedHmmPane.setTitleAt(0, hmm.name);
		hmmText.setHmm(hmm);
		hmm.Plan7Renormalize();
		hmmString = hmm.consensus;
		hmmLayout = new HmmLayout(hmm, 20);
		hmmCanvas = new HmmPanel(hmmLayout, hmm, jChartLabel,
				commandHistory, jHmmScrollPane);
		hmmCanvas.setBackground(Color.WHITE);

		nullLayout = new NullModelLayout(hmm, 20);
		nullCanvas = new NullModelPanel(nullLayout, hmm,
				jChartLabel, commandHistory, jNullScrollPane);
		nullCanvas.setBackground(Color.WHITE);

		logoCanvas = new LogoPanel(hmm, 20, jLogoScrollPane);
		logoCanvas.setBackground(Color.WHITE);

		jHmmScrollPane.setViewportView(hmmCanvas);
		jNullScrollPane.setViewportView(nullCanvas);
		jLogoScrollPane.setViewportView(logoCanvas);

		jChartLabel.setAlphabet(hmm.alphabets);
		commandHistory.clear();
		tabbedHmmPane.setSelectedIndex(0);
	}
	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("OpenHmm")) {
			if (hmm != null && hmm.modified) {
				int ret = JOptionPane.showConfirmDialog(this,
						"Hmm model has been modified, do you want to save?",
						"Confirm", JOptionPane.YES_NO_CANCEL_OPTION);
				if (ret == JOptionPane.YES_OPTION) {
					ActionEvent e1 = new ActionEvent(this,
							ActionEvent.ACTION_PERFORMED, "SaveHmm");
					actionPerformed(e1);
				} else if (ret == JOptionPane.CANCEL_OPTION)
					return;
			}
			JFileChooser fc = new JFileChooser(currentDirectory);
			fc.addChoosableFileFilter(new HmmFileFilter());
			int returnVal = fc.showDialog(this, "Open");
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				hmmFileName = fc.getSelectedFile().getAbsolutePath();
				hmm = new Hmm();
				try {
					hmm.read(hmmFileName);
					initWindowLoadHmm();
				} catch (IOException e1) {
					JOptionPane.showMessageDialog(this,
							"Can not read hmm model file");
				} catch (InvalidFormatException e1) {
					JOptionPane.showMessageDialog(this, "File format error");
				}
			}
            currentDirectory=fc.getCurrentDirectory();

		}
		if (e.getActionCommand().equals("OpenPfam")) {
			if (hmm != null && hmm.modified) {
				int ret = JOptionPane.showConfirmDialog(this,
						"Hmm model has been modified, do you want to save?",
						"Confirm", JOptionPane.YES_NO_CANCEL_OPTION);
				if (ret == JOptionPane.YES_OPTION) {
					ActionEvent e1 = new ActionEvent(this,
							ActionEvent.ACTION_PERFORMED, "SaveHmm");
					actionPerformed(e1);
				} else if (ret == JOptionPane.CANCEL_OPTION)
					return;
			}
			JFileChooser fc = new JFileChooser(currentDirectory);
			int returnVal = fc.showDialog(this, "Open Pfam");
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				pfamFileName = fc.getSelectedFile().getAbsolutePath();
				
				if (pfam==null || !pfamFileName.equals(pfam.pfamFileName))
				{
					try
					{
						pfam = new Pfam(pfamFileName);
					}
					catch (IOException e1)
					{
						JOptionPane.showMessageDialog(this, "Can not read pfam file");
						return;
					}
				}
				PfamPickDialog pfamPanel = new PfamPickDialog(pfam);
				pfamPanel.setVisible(true);
				if (pfamPanel.selectedHmmIndex!=-1)
				{
					try {
						hmm=pfam.getHmm(pfamPanel.selectedHmmIndex);
						initWindowLoadHmm();
					} catch (IOException e1) {
						JOptionPane.showMessageDialog(this,
								"Can not read hmm model file");
					} catch (InvalidFormatException e1) {
						JOptionPane.showMessageDialog(this, "File format error");
					}
				}
			}
            currentDirectory=fc.getCurrentDirectory();
		}
		if (e.getActionCommand().equals("SaveHmm")) {
			if (hmm == null) {
				JOptionPane.showMessageDialog(this, "Nothing to save");
				return;
			}
			JFileChooser fc = new JFileChooser(currentDirectory);
			fc.addChoosableFileFilter(new HmmFileFilter());
			int returnVal = fc.showDialog(this, "Save");
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				hmmFileName = fc.getSelectedFile().getAbsolutePath();
				if (!hmmFileName.toUpperCase().endsWith(".HMM"))
					hmmFileName += ".hmm";
				if (new File(hmmFileName).exists()) {
					int ret = JOptionPane.showConfirmDialog(this,
							"Do you want to overwrite existing file?");
					if (ret != JOptionPane.YES_OPTION) {
						return;
					}
				}
				try {
					hmm.write(hmmFileName);
				} catch (IOException e1) {
					JOptionPane.showMessageDialog(this,
							"Can not write hmm model file");
				}
				hmm.modified = false;
			}
            currentDirectory=fc.getCurrentDirectory();
		} else if (e.getActionCommand().equals("OpenSequence")) {
			JFileChooser fc = new JFileChooser(currentDirectory);
			fc.addChoosableFileFilter(new SequenceFileFilter());
			int returnVal = fc.showDialog(this, "Open");
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				sequenceFileName = fc.getSelectedFile().getAbsolutePath();
				SQFILE sqfile = new SQFILE();
				try {
					sequences = sqfile.ReadMultipleRseqs(sequenceFileName,
							SQFILE.SQFILE_FORMAT.SQFILE_UNKNOWN);
					Object[][] data = new Object[sequences.size()][2];
					for (int i = 0; i < sequences.size(); i++) {
						data[i][0] = sequences.get(i).sqInfo.name;
						data[i][1] = sequences.get(i).seqs;
					}
					Object[] columns = new Object[] { "Name", "Sequence" };
					TableModel model = new DefaultTableModel(data, columns) {
						public boolean isCellEditable(int row, int col) {
							return false;
						}
					};
					jSequenceTable.setModel(model);
					jSequenceTable.getColumnModel().getColumn(0)
							.setPreferredWidth(100);
					jSequenceTable.getColumnModel().getColumn(1)
							.setPreferredWidth(500);
				} catch (Exception e1) {
					JOptionPane.showMessageDialog(this,
							"Can not read sequence file");
				}
			}
            currentDirectory=fc.getCurrentDirectory();
		} else if (e.getActionCommand().equals("SaveImage")) {
			class JpgFileFilter extends FileFilter {
				public boolean accept(File f) {
					if (f.isDirectory()) {
						return true;
					} else if (f.getName().toUpperCase().endsWith("JPG"))
						return true;
					return false;
				}

				public String getDescription() {
					return "JPG File(*.jpg)";
				}
			}
			;
			class PngFileFilter extends FileFilter {
				public boolean accept(File f) {
					if (f.isDirectory()) {
						return true;
					} else if (f.getName().toUpperCase().endsWith("PNG"))
						return true;
					return false;
				}

				public String getDescription() {
					return "PNG File(*.png)";
				}
			}
			;
			if (tabbedHmmPane.getSelectedIndex() == 0 && hmmCanvas == null
					|| tabbedHmmPane.getSelectedIndex() == 1
					&& nullCanvas == null
					|| tabbedHmmPane.getSelectedIndex() == 2
					&& logoCanvas == null
					|| tabbedHmmPane.getSelectedIndex() == 3) {
				JOptionPane.showMessageDialog(this, "Nothing to save");
				return;
			}
			JFileChooser fc = new JFileChooser(currentDirectory);
			fc.addChoosableFileFilter(new PngFileFilter());
			fc.addChoosableFileFilter(new JpgFileFilter());
			int returnVal = fc.showDialog(this, "Save");
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				String imageFileName = fc.getSelectedFile().getAbsolutePath();
				String imageType;
				if (fc.getFileFilter() instanceof JpgFileFilter) {
					imageType = "jpg";
					if (!imageFileName.toUpperCase().endsWith(".JPG")) {
						imageFileName += ".jpg";
					}
				} else if (fc.getFileFilter() instanceof PngFileFilter) {
					imageType = "png";
					if (!imageFileName.toUpperCase().endsWith(".PNG")) {
						imageFileName += ".png";
					}
				} else {
					if (!imageFileName.toUpperCase().endsWith(".PNG")) {
						imageType = "png";
					} else if (!imageFileName.toUpperCase().endsWith(".JPG")) {
						imageType = "jpg";
					} else {
						JOptionPane.showMessageDialog(this,
								"Unsupported image file format");
						return;
					}
				}
				if (new File(imageFileName).exists()) {
					int ret = JOptionPane.showConfirmDialog(this,
							"Do you want to overwrite existing file?");
					if (ret != JOptionPane.YES_OPTION) {
						return;
					}
				}
				BufferedImage image = null;
				if (tabbedHmmPane.getSelectedIndex() == 0) {
					int width = hmmCanvas.getWidth();
					int height = hmmCanvas.getHeight();
					image = new BufferedImage(width, height,
							BufferedImage.TYPE_INT_RGB);
					Graphics2D g2 = image.createGraphics();
					hmmCanvas.bitmapMode = true;
					hmmCanvas.paint(g2);
					hmmCanvas.bitmapMode = false;
					g2.dispose();
				} else if (tabbedHmmPane.getSelectedIndex() == 1) {
					int width = nullCanvas.getWidth();
					int height = nullCanvas.getHeight();
					image = new BufferedImage(width, height,
							BufferedImage.TYPE_INT_RGB);
					Graphics2D g2 = image.createGraphics();
					nullCanvas.paint(g2);
					g2.dispose();
				} else if (tabbedHmmPane.getSelectedIndex() == 2) {
					int width = logoCanvas.getWidth();
					int height = logoCanvas.getHeight();
					image = new BufferedImage(width, height,
							BufferedImage.TYPE_INT_RGB);
					Graphics2D g2 = image.createGraphics();
					logoCanvas.bitmapMode = true;
					logoCanvas.paint(g2);
					logoCanvas.bitmapMode = false;
					g2.dispose();
				}
				try {
					ImageIO.write(image, imageType, new File(imageFileName));
				} catch (IOException ioe) {
					System.out.println(ioe.getMessage());
				}
			}
            currentDirectory=fc.getCurrentDirectory();
		} else if (e.getActionCommand().equals("Undo")) {
			HmmCommand.Target target = commandHistory.undo();
			if (target == HmmCommand.Target.ProfileHmm
					|| target == HmmCommand.Target.NNode
					|| target == HmmCommand.Target.BNode
					|| target == HmmCommand.Target.ENode
					|| target == HmmCommand.Target.CNode
					|| target == HmmCommand.Target.JNode) {
				hmmLayout.setHmmValue(hmm);
				hmmLayout.generateShape(hmmLayout.size);
				hmmCanvas.repaint();
			} else {
				nullLayout.setHmmValue(hmm);
				nullLayout.generateShape(nullLayout.size);
				nullCanvas.repaint();
			}
		} else if (e.getActionCommand().equals("Redo")) {
			HmmCommand.Target target = commandHistory.redo();
			if (target == HmmCommand.Target.ProfileHmm
					|| target == HmmCommand.Target.NNode
					|| target == HmmCommand.Target.BNode
					|| target == HmmCommand.Target.ENode
					|| target == HmmCommand.Target.CNode
					|| target == HmmCommand.Target.JNode) {
				hmmLayout.setHmmValue(hmm);
				hmmLayout.generateShape(hmmLayout.size);
				hmmCanvas.repaint();
			} else {
				nullLayout.setHmmValue(hmm);
				nullLayout.generateShape(nullLayout.size);
				nullCanvas.repaint();
			}
		} else if (e.getActionCommand().equals("AlignPath")) {
			if (hmm == null) {
				JOptionPane.showMessageDialog(this, "Please open one hmm file");
				return;
			}
			if (sequences == null) {
				JOptionPane.showMessageDialog(this,
						"Please open one sequence file");
				return;
			}
			if (jSequenceTable.getSelectedRow() == -1) {
				JOptionPane.showMessageDialog(this,
						"Please select one sequence");
				return;
			}
			if(hmm.hmmerVersion==3) {
				JOptionPane.showMessageDialog(this, "Cannot align to HMMER3 models");
				return;
			}
			String seqs = (String) jSequenceTable.getModel().getValueAt(
					jSequenceTable.getSelectedRow(), 1);
			P7Trace tr = hmm.doAlign(seqs);
			hmmLayout.tagPath(tr, seqs);
			hmmCanvas.repaint();
		} else if (e.getActionCommand().equals("AlignAll")) {
			if (hmm == null) {
				JOptionPane.showMessageDialog(this, "Please open one hmm file");
				return;
			}
			if (sequences == null) {
				JOptionPane.showMessageDialog(this,
						"Please open one sequence file");
				return;
			}
			if(hmm.hmmerVersion==3) {
				JOptionPane.showMessageDialog(this, "Cannot align to HMMER3 models");
				return;
			}
			msa = hmm.alignAll(sequences);
			Object[][] data = new Object[sequences.size()][3];
			for (int i = 0; i < sequences.size(); i++) {
				data[i][0] = msa.sqname.get(i);
				data[i][1] = msa.scores[i];
				data[i][2] = msa.aseq.get(i);
			}
			Object[] columns = new Object[] { "Name", "Score", "Sequence" };
			TableModel model = new DefaultTableModel(data, columns) {
				public boolean isCellEditable(int row, int column) {
					return false;
				}
			};
			alignResultTable.setModel(model);
			alignResultTable.getColumnModel().getColumn(0).setPreferredWidth(
					100);
			alignResultTable.getColumnModel().getColumn(1).setPreferredWidth(
					100);
			alignResultTable.getColumnModel().getColumn(2).setPreferredWidth(
					500);
			tabbedInfoPane.setSelectedIndex(1);
		} else if (e.getActionCommand().equals("SaveSequence")) {
			if (msa == null) {
				JOptionPane.showMessageDialog(this, "No Sequence to save");
				return;
			}

			if (msa == null) {
				JOptionPane.showMessageDialog(this, "Nothing to save");
				return;
			}
			JFileChooser fc = new JFileChooser(currentDirectory);
			fc.addChoosableFileFilter(new MSFFileFilter());
			int returnVal = fc.showDialog(this, "Save");
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				sequenceFileName = fc.getSelectedFile().getAbsolutePath();
				if (new File(sequenceFileName).exists()) {
					int ret = JOptionPane.showConfirmDialog(this,
							"Do you want to overwrite existing file?");
					if (ret != JOptionPane.YES_OPTION) {
						return;
					}
				}
				try {
					MSAFILE msaFile = new MSAFILE();
					FileWriter w = new FileWriter(sequenceFileName);
					msaFile.write(w, msa, SQFILE.SQFILE_FORMAT.MSAFILE_MSF,
							false);
				} catch (IOException e1) {
					JOptionPane.showMessageDialog(this,
							"Can not write sequence file");
				}
			}
            currentDirectory=fc.getCurrentDirectory();
		} else if (e.getActionCommand().equals("Clear")) {
			if (hmm != null) {
				hmmLayout.clearPath();
				hmmCanvas.repaint();
			}
		} else if (e.getActionCommand().equals("ZoomIn")) {
			if (tabbedHmmPane.getSelectedIndex() == 0) {
				if (hmmLayout != null) {
					hmmLayout.resize(ZoomFactor);
					if (hmmCanvas != null) {
						hmmCanvas.resize();
						hmmCanvas.repaint();
						hmmCanvas.revalidate();
					}
				}
			} else if (tabbedHmmPane.getSelectedIndex() == 1) {
				if (nullLayout != null) {
					nullLayout.resize(ZoomFactor);
					if (nullCanvas != null) {
						nullCanvas.resize();
						nullCanvas.repaint();
						nullCanvas.revalidate();
					}
				}
			} else if (tabbedHmmPane.getSelectedIndex() == 2) {
				if (logoCanvas != null) {
					logoCanvas.resize(ZoomFactor);
					logoCanvas.repaint();
					logoCanvas.revalidate();
				}
			}
		} else if (e.getActionCommand().equals("ZoomOut")) {
			if (tabbedHmmPane.getSelectedIndex() == 0) {
				if (hmmLayout != null) {
					hmmLayout.resize((double) 1 / ZoomFactor);
					if (hmmCanvas != null) {
						hmmCanvas.resize();
						hmmCanvas.repaint();
						hmmCanvas.revalidate();
					}
				}
			} else if (tabbedHmmPane.getSelectedIndex() == 1) {
				if (nullLayout != null) {
					nullLayout.resize((double) 1 / ZoomFactor);
					if (nullCanvas != null) {
						nullCanvas.resize();
						nullCanvas.repaint();
						nullCanvas.revalidate();
					}
				}
			} else if (tabbedHmmPane.getSelectedIndex() == 2) {
				if (logoCanvas != null) {
					logoCanvas.resize((double) 1 / ZoomFactor);
					logoCanvas.repaint();
					logoCanvas.revalidate();
				}
			}

		} else if (e.getActionCommand().equals("IncreasePrecision")) {
			if (hmm != null && precision < MaxPrecision) {
				precision++;
				hmmLayout.precision = precision;
				nullLayout.precision = precision;
				hmmCanvas.repaint();
				nullCanvas.repaint();
			}
		} else if (e.getActionCommand().equals("DecreasePrecision")) {
			if (hmm != null && precision > MinPrecision) {
				precision--;
				hmmLayout.precision = precision;
				nullLayout.precision = precision;
				hmmCanvas.repaint();
				nullCanvas.repaint();
			}
		} else if (e.getActionCommand().equals("DefaultLayout")) {
			if (tabbedHmmPane.getSelectedIndex() == 0) {
				if (hmmLayout != null) {
					hmmLayout.defaultLayout();
					hmmLayout.generateShape(hmmLayout.size);
					hmmCanvas.repaint();
					hmmCanvas.revalidate();
				}
			} else if (tabbedHmmPane.getSelectedIndex() == 1) {
				if (nullLayout != null) {
					nullLayout.defaultLayout();
					nullLayout.generateShape(nullLayout.size);
					nullCanvas.repaint();
					nullCanvas.revalidate();
				}
			}
		} else if (e.getActionCommand().equals("ReNumber")) {
			if (tabbedHmmPane.getSelectedIndex() == 0) {
				if (hmmLayout != null) {
					hmmLayout.reNumber();
					hmmCanvas.repaint();
					hmmCanvas.revalidate();
				}
			}
		} else if (e.getActionCommand().equals("ViewLabel")) {
			if (hmmLayout != null) {
				for (int i = 1; i <= hmmLayout.M; i++) {
					hmmLayout.mNodes[i].showConsensus = false;
				}
				hmmCanvas.repaint();
			}
		} else if (e.getActionCommand().equals("ViewConsensus")) {
			if (hmmLayout != null) {
				for (int i = 1; i <= hmmLayout.M; i++) {
					hmmLayout.mNodes[i].showConsensus = true;
				}
				hmmCanvas.repaint();
			}
		} else if (e.getActionCommand().equals("HidePath")) {
			if (hmm != null) {
				if (hidePathMenuItem.getState() == true) {
					for (int i = 1; i <= hmm.M; i++) {
						if (i != 1 && hmm.begin[i] == 0)
							hmmLayout.bm[i].show = false;
						if (i != hmm.M && hmm.end[i] == 0)
							hmmLayout.me[i].show = false;
					}
					if (hmm.xt[Hmm.XTN][Hmm.LOOP] == 0)
						hmmLayout.nn.show = false;
					if (hmm.xt[Hmm.XTC][Hmm.LOOP] == 0)
						hmmLayout.cc.show = false;
					if (hmm.xt[Hmm.XTE][Hmm.LOOP] == 0) {
						hmmLayout.jNode.show = false;
						hmmLayout.ej.show = false;
						hmmLayout.jb.show = false;
						hmmLayout.jj.show = false;
					}
				} else {
					for (int i = 1; i <= hmm.M; i++) {
						hmmLayout.me[i].show = true;
						hmmLayout.bm[i].show = true;
					}
					hmmLayout.nn.show = true;
					hmmLayout.cc.show = true;
					hmmLayout.jNode.show = true;
					hmmLayout.ej.show = true;
					hmmLayout.jb.show = true;
					hmmLayout.jj.show = true;
				}
				hmmCanvas.repaint();
			}
		} else if (e.getActionCommand().equals("HideEmission")) {
			if (hideEmissionMenuItem.getState() == true)
				hmmLayout.showEmission = false;
			else
				hmmLayout.showEmission = true;
			hmmCanvas.repaint();
		} else if (e.getActionCommand().equals("Exit")) {
			if (hmm == null || !hmm.modified)
				System.exit(0);
			int ret = JOptionPane.showConfirmDialog(this,
					"Hmm model has been modified, do you want to save?",
					"Confirm", JOptionPane.YES_NO_CANCEL_OPTION);
			if (ret == JOptionPane.YES_OPTION) {
				ActionEvent e1 = new ActionEvent(this,
						ActionEvent.ACTION_PERFORMED, "SaveHmm");
				actionPerformed(e1);
				System.exit(0);
			} else if (ret == JOptionPane.NO_OPTION) {
				System.exit(0);
			}
		}
	}

	public void componentHidden(ComponentEvent e) {
	}

	public void componentMoved(ComponentEvent e) {
	}

	public void componentResized(ComponentEvent e) {
		resize();
	}

	public void componentShown(ComponentEvent e) {
	}

	public void windowActivated(WindowEvent e) {
	}

	public void windowClosed(WindowEvent e) {
	}

	public void windowClosing(WindowEvent e) {
		ActionEvent e1 = new ActionEvent(this, ActionEvent.ACTION_PERFORMED,
				"Exit");
		actionPerformed(e1);
	}

	public void windowDeactivated(WindowEvent e) {
	}

	public void windowDeiconified(WindowEvent e) {
	}

	public void windowIconified(WindowEvent e) {
	}

	public void windowOpened(WindowEvent e) {
	}

	public void mouseClicked(MouseEvent e) {
		if (e.getSource() == jSequenceTable) {
			if (e.getClickCount() == 2 && jSequenceTable.getSelectedRow() != -1
					&& hmm != null) {
				SequenceEditor editor = new SequenceEditor(
						(String) jSequenceTable.getModel().getValueAt(
								jSequenceTable.getSelectedRow(), 1),
						hmm.alphabet, true);
				editor.setVisible(true);
				if (editor.text != null)
					jSequenceTable.getModel().setValueAt(editor.text,
							jSequenceTable.getSelectedRow(), 1);
			}
		} else if (e.getSource() == alignResultTable) {
			if (e.getClickCount() == 2
					&& alignResultTable.getSelectedRow() != -1 && hmm != null) {
				SequenceEditor editor = new SequenceEditor(
						((StringBuffer) alignResultTable.getModel().getValueAt(
								alignResultTable.getSelectedRow(), 2))
								.toString(), null, false);
				editor.setVisible(true);
			}
		}
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	public void mousePressed(MouseEvent e) {
	}

	public void mouseReleased(MouseEvent e) {
	}

	public void stateChanged(ChangeEvent e) {
		if (e.getSource() == tabbedHmmPane) {
			int sel = tabbedHmmPane.getSelectedIndex();
			if (sel == 2 && hmm != null) {
				logoCanvas.generateLogo();
				logoCanvas.repaint();
			}
			if (sel == 3 && hmm != null) {
				hmmText.updateHmmText();
				hmmText.setCaretPosition(0);
			}
		}
	}
}
