package edu.ucf.cs.hmm.visualizer;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Rectangle2D;

import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;

/**
 * Panel to show Hmm layout
 * 
 */
public class HmmPanel extends JPanel implements MouseMotionListener,
		MouseListener, ActionListener {
	Hmm hmm;

	HmmLayout layout;

	NodeState lingerShape = null;

	NodeState selectedState = null;

	NodeState actionTarget = null;

	Point selectedPosition = null;

	Point lastDragDrawPosition = null;

	StateChartLabel stateChartLabel;

	JPopupMenu popupMenu1, popupMenu2;

	JMenuItem dupBeforeMenuItem;

	JMenuItem dupAfterMenuItem;

	HmmCommandList commandHistory;

	Rectangle2D.Float clip;

	JScrollPane container;

	boolean bitmapMode = false;

	public HmmPanel(HmmLayout _layout, Hmm _hmm,
			StateChartLabel _stateChartLabel, HmmCommandList _commandHistory,
			JScrollPane _container) {
		layout = _layout;
		hmm = _hmm;
		commandHistory = _commandHistory;
		stateChartLabel = _stateChartLabel;
		container = _container;
		setPreferredSize(layout.getDimension());
		addMouseMotionListener(this);
		addMouseListener(this);
		setAutoscrolls(true);
		popupMenu1 = new JPopupMenu();

		JMenuItem menuItem = new JMenuItem("Modify Node...", KeyEvent.VK_E);
		menuItem.setActionCommand("Modify");
		menuItem.addActionListener(this);
		popupMenu1.add(menuItem);

		menuItem = new JMenuItem("Delete...", KeyEvent.VK_D);
		menuItem.setActionCommand("Delete");
		menuItem.addActionListener(this);
		popupMenu1.add(menuItem);

		dupBeforeMenuItem = new JMenuItem("Duplicate State Before",
				KeyEvent.VK_B);
		dupBeforeMenuItem.setActionCommand("Before");
		dupBeforeMenuItem.addActionListener(this);
		popupMenu1.add(dupBeforeMenuItem);

		dupAfterMenuItem = new JMenuItem("Duplicate State After", KeyEvent.VK_E);
		dupAfterMenuItem.setActionCommand("After");
		dupAfterMenuItem.addActionListener(this);
		popupMenu1.add(dupAfterMenuItem);

		popupMenu2 = new JPopupMenu();
		menuItem = new JMenuItem("Modify Node...", KeyEvent.VK_E);
		menuItem.setActionCommand("Modify");
		menuItem.addActionListener(this);
		popupMenu2.add(menuItem);
	}

	public void resize() {
		setPreferredSize(layout.getDimension());
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Graphics2D g2d = (Graphics2D) g;
		Rectangle2D.Float clip;
		if (!bitmapMode) {
			clip = new Rectangle2D.Float();
			clip.x = container.getViewport().getViewPosition().x;
			clip.y = container.getViewport().getViewPosition().y;
			clip.width = container.getViewport().getExtentSize().width;
			clip.height = container.getViewport().getExtentSize().height;
		} else
			clip = null;
		layout.drawWithLabel(g2d, clip);
	}

	/**
	 * Show the shape of the node when mouse dragging
	 */
	public void mouseDragged(MouseEvent e) {
		Rectangle r = new Rectangle(e.getX(), e.getY(), 1, 1);
		scrollRectToVisible(r);
		if (selectedState != null) {
			int offsetx = e.getPoint().x - selectedState.position.x;
			int offsety = e.getPoint().y - selectedState.position.y;

			Graphics2D g2d = (Graphics2D) getGraphics();
			g2d.setXORMode(Color.GRAY);
			if (lastDragDrawPosition != null) {
				Shape l = selectedState.generateShapeEx(lastDragDrawPosition.x,
						lastDragDrawPosition.y);
				g2d.draw(l);
			}
			Shape s = selectedState.generateShapeEx(selectedState.position.x
					+ offsetx, selectedState.position.y + offsety);
			g2d.draw(s);
			lastDragDrawPosition = new Point(
					selectedState.position.x + offsetx,
					selectedState.position.y + offsety);
		}
	}

	/**
	 * Draw emission probability when mouse is inside applicable node
	 */
	public void mouseMoved(MouseEvent e) {
		if (popupMenu1.isVisible() || popupMenu2.isVisible())
			return;
		if (hmm != null) {
			Graphics2D g2d = (Graphics2D) getGraphics();
			if (lingerShape != null && lingerShape.contains(e.getPoint()))
				return;
			if (lingerShape != null) {
				lingerShape.toggleSelected(g2d, false);
				lingerShape = null;
			}
			for (int i = 1; i <= hmm.M; i++) {
				if (layout.mNodes[i].contains(e.getPoint())) {
					lingerShape = layout.mNodes[i];
					lingerShape.toggleSelected(g2d, true);
					stateChartLabel.setState("M" + i, hmm.mat[i]);
					stateChartLabel.repaint();
					return;
				} else if (i != hmm.M
						&& layout.iNodes[i].contains(e.getPoint())) {
					lingerShape = layout.iNodes[i];
					lingerShape.toggleSelected(g2d, true);
					stateChartLabel.setState("I" + i, hmm.ins[i]);
					stateChartLabel.repaint();
					return;
				}
			}
			if (layout.nn.contains(e.getPoint())) {
				stateChartLabel.setState("NN", hmm.nullModel);
				stateChartLabel.repaint();
			} else if (layout.cc.contains(e.getPoint())) {
				stateChartLabel.setState("CC", hmm.nullModel);
				stateChartLabel.repaint();
			} else if (layout.jj.contains(e.getPoint())) {
				stateChartLabel.setState("JJ", hmm.nullModel);
				stateChartLabel.repaint();
			}
		}
	}

	/**
	 * Double click or right click
	 */
	public void mouseClicked(MouseEvent e) {
		// Double click
		if (e.getButton() == MouseEvent.BUTTON1 && e.getClickCount() != 2
				|| e.getButton() == MouseEvent.BUTTON3
				&& e.getClickCount() != 1) {
			return;
		}
		ActionEvent ae = new ActionEvent(popupMenu1,
				ActionEvent.ACTION_PERFORMED, "Modify");
		for (int i = 1; i <= hmm.M; i++) {
			if (layout.mNodes[i].contains(e.getPoint())) {
				actionTarget = layout.mNodes[i];
				if (e.getButton() == MouseEvent.BUTTON3)
					popupMenu1.show(e.getComponent(), e.getX(), e.getY());
				else
					actionPerformed(ae);
				return;
			} else if (layout.dNodes[i].contains(e.getPoint())) {
				actionTarget = layout.dNodes[i];
				if (e.getButton() == MouseEvent.BUTTON3)
					popupMenu2.show(e.getComponent(), e.getX(), e.getY());
				else
					actionPerformed(ae);
				return;
			} else if (i != hmm.M && layout.iNodes[i].contains(e.getPoint())) {
				actionTarget = layout.iNodes[i];
				if (e.getButton() == MouseEvent.BUTTON3)
					popupMenu2.show(e.getComponent(), e.getX(), e.getY());
				else
					actionPerformed(ae);
				return;
			}
		}
		if (layout.nNode.contains(e.getPoint())) {
			actionTarget = layout.nNode;
			if (e.getButton() == MouseEvent.BUTTON3)
				popupMenu2.show(e.getComponent(), e.getX(), e.getY());
			else
				actionPerformed(ae);
		} else if (layout.bNode.contains(e.getPoint())) {
			actionTarget = layout.bNode;
			if (e.getButton() == MouseEvent.BUTTON3)
				popupMenu2.show(e.getComponent(), e.getX(), e.getY());
			else
				actionPerformed(ae);
		} else if (layout.eNode.contains(e.getPoint())) {
			actionTarget = layout.eNode;
			if (e.getButton() == MouseEvent.BUTTON3)
				popupMenu2.show(e.getComponent(), e.getX(), e.getY());
			else
				actionPerformed(ae);

		} else if (layout.cNode.contains(e.getPoint())) {
			actionTarget = layout.cNode;
			if (e.getButton() == MouseEvent.BUTTON3)
				popupMenu2.show(e.getComponent(), e.getX(), e.getY());
			else
				actionPerformed(ae);

		} else if (layout.jNode.contains(e.getPoint())) {
			actionTarget = layout.jNode;
			if (e.getButton() == MouseEvent.BUTTON3)
				popupMenu2.show(e.getComponent(), e.getX(), e.getY());
			else
				actionPerformed(ae);
		}
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	/**
	 * Prepare for dragging if applicable
	 */
	public void mousePressed(MouseEvent e) {
		if (e.getButton() != MouseEvent.BUTTON1)
			return;
		NodeState state = null;
		for (int i = 1; i <= hmm.M; i++) {
			if (layout.mNodes[i].contains(e.getPoint())) {
				state = layout.mNodes[i];
				break;
			} else if (i != hmm.M && layout.iNodes[i].contains(e.getPoint())) {
				state = layout.iNodes[i];
				break;
			} else if (layout.dNodes[i].contains(e.getPoint())) {
				state = layout.dNodes[i];
				break;
			}
		}
		if (state == null) {
			if (layout.sNode.contains(e.getPoint())) {
				state = layout.sNode;
			} else if (layout.nNode.contains(e.getPoint())) {
				state = layout.nNode;
			} else if (layout.bNode.contains(e.getPoint())) {
				state = layout.bNode;
			} else if (layout.eNode.contains(e.getPoint())) {
				state = layout.eNode;
			} else if (layout.cNode.contains(e.getPoint())) {
				state = layout.cNode;
			} else if (layout.tNode.contains(e.getPoint())) {
				state = layout.tNode;
			} else if (layout.jNode.contains(e.getPoint())) {
				state = layout.jNode;
			}
		}
		if (state != null) {
			selectedState = state;
			selectedPosition = e.getPoint();
		} else {
			clearSelection();
		}
	}

	/**
	 * Clear previous dragging node
	 *
	 */
	public void clearSelection() {
		selectedState = null;
		selectedPosition = null;
		lastDragDrawPosition = null;
		lingerShape = null;
	}

	/**
	 * Release dragging
	 */
	public void mouseReleased(MouseEvent e) {
		if (e.getButton() != MouseEvent.BUTTON1)
			return;
		if (lastDragDrawPosition != null) {
			layout.moveNode(selectedState, lastDragDrawPosition.x,
					lastDragDrawPosition.y);
			repaint();
		}
		selectedState = null;
		selectedPosition = null;
		lastDragDrawPosition = null;
	}

	/**
	 * Command handler
	 */
	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("Modify")) {
			Probabilities tp;
			Probabilities ep = new Probabilities(hmm.alphabets);
			int position = -1;
			boolean modified = false;
			for (int i = 1; i <= hmm.M; i++) {
				if (actionTarget == layout.mNodes[i]
						|| (i != hmm.M && actionTarget == layout.iNodes[i])
						|| actionTarget == layout.dNodes[i]) {
					position = i;
				}
			}
			if (position != -1) // something selected
			{
				HmmCommand command;
				command = new HmmCommand(hmm, layout, null);
				command.target = HmmCommand.Target.ProfileHmm;
				command.commandType = HmmCommand.CommandType.Modify;
				command.position = position;
				command.oldState.mat = hmm.mat[position].clone();
				command.oldState.begin = hmm.begin[position];
				command.oldState.end = hmm.end[position];
				command.oldState.modified[0] = layout.mNodes[position].modified;
				command.oldState.modified[2] = layout.dNodes[position].modified;
				if (position != hmm.M) {
					command.oldState.modified[1] = layout.iNodes[position].modified;
					command.oldState.t = hmm.t[position].clone();
					command.oldState.ins = hmm.ins[position].clone();
				}

				if (actionTarget == layout.mNodes[position]) // Modify m node
				{
					if (position != hmm.M) {
						tp = new Probabilities(4);
						tp.add("M" + (position + 1), hmm.t[position][0]);
						tp.add("I" + position, hmm.t[position][1]);
						tp.add("D" + (position + 1), hmm.t[position][2]);
						tp.add("E", hmm.end[position]);
					} else {
						tp = new Probabilities(1);
						tp.add("E", 1.0);
					}
					ep.set(hmm.mat[position]);
					JModifyNodeDialog dialog = new JModifyNodeDialog(
							actionTarget.label, tp, ep);
					dialog.setVisible(true);
					if (dialog.action == JOptionPane.OK_OPTION) {
						modified = true;
						hmm.mat[position] = ep.probs;
						actionTarget.modified = true;
						hmm.modified = true;
						if (position != hmm.M) {
							hmm.t[position][0] = tp.probs[0];
							hmm.t[position][1] = tp.probs[1];
							hmm.t[position][2] = tp.probs[2];
							hmm.end[position] = tp.probs[3];
						} else
							hmm.end[position] = tp.probs[0];
						hmm.getRepresentativeString();
						layout.setHmmValue(hmm);
						layout.generateShape(layout.size);
						repaint();
					}
				} else if (position != hmm.M
						&& actionTarget == layout.iNodes[position]) // Modify i
																	// node
				{
					tp = new Probabilities(2);
					tp.add("M" + (position + 1), hmm.t[position][3]);
					tp.add("I" + position, hmm.t[position][4]);
					ep.set(hmm.ins[position]);
					JModifyNodeDialog dialog = new JModifyNodeDialog(
							actionTarget.label, tp, ep);
					dialog.setVisible(true);
					if (dialog.action == JOptionPane.OK_OPTION) {
						modified = true;
						hmm.ins[position] = ep.probs;
						actionTarget.modified = true;
						hmm.modified = true;
						hmm.t[position][3] = tp.probs[0];
						hmm.t[position][4] = tp.probs[1];
						layout.setHmmValue(hmm);
						layout.generateShape(layout.size);
						repaint();
					}
					layout.generateShape(layout.size);
				} else if (actionTarget == layout.dNodes[position]) // Modify d
																	// node
				{
					if (position != hmm.M) {
						tp = new Probabilities(2);
						tp.add("M" + (position + 1), hmm.t[position][5]);
						tp.add("D" + (position + 1), hmm.t[position][6]);

						JModifyNodeDialog dialog = new JModifyNodeDialog(
								actionTarget.label, tp, null);
						dialog.setVisible(true);
						if (dialog.action == JOptionPane.OK_OPTION) {
							modified = true;
							actionTarget.modified = true;
							hmm.modified = true;
							if (position != hmm.M) {
								hmm.t[position][5] = tp.probs[0];
								hmm.t[position][6] = tp.probs[1];
							}
							layout.setHmmValue(hmm);
							layout.generateShape(layout.size);
							repaint();
						}
					}
				}
				// save node to command history
				if (modified) {
					command.newState.mat = hmm.mat[position].clone();
					command.newState.begin = hmm.begin[position];
					command.newState.end = hmm.end[position];
					command.newState.modified[0] = layout.mNodes[position].modified;
					command.newState.modified[2] = layout.dNodes[position].modified;
					if (position != hmm.M) {
						command.newState.modified[1] = layout.iNodes[position].modified;
						command.newState.t = hmm.t[position].clone();
						command.newState.ins = hmm.ins[position];
					}
					commandHistory.add(command);
				}
			}
			if (actionTarget == layout.nNode) {
				tp = new Probabilities(2);
				tp.add("B", hmm.xt[Hmm.XTN][Hmm.MOVE]);
				tp.add("N", hmm.xt[Hmm.XTN][Hmm.LOOP]);
				JModifyNodeDialog dialog = new JModifyNodeDialog("N", tp, null);
				dialog.setVisible(true);
				if (dialog.action == JOptionPane.OK_OPTION) {
					HmmCommand command = new HmmCommand(hmm, layout, null);
					command.target = HmmCommand.Target.NNode;
					command.commandType = HmmCommand.CommandType.Modify;
					command.oldState.xt[Hmm.XTN][Hmm.MOVE] = hmm.xt[Hmm.XTN][Hmm.MOVE];
					command.oldState.xt[Hmm.XTN][Hmm.LOOP] = hmm.xt[Hmm.XTN][Hmm.LOOP];
					command.oldState.modified[0] = actionTarget.modified;
					actionTarget.modified = true;
					hmm.modified = true;
					hmm.xt[Hmm.XTN][Hmm.MOVE] = tp.probs[0];
					hmm.xt[Hmm.XTN][Hmm.LOOP] = tp.probs[1];
					command.newState.xt[Hmm.XTN][Hmm.MOVE] = hmm.xt[Hmm.XTN][Hmm.MOVE];
					command.newState.xt[Hmm.XTN][Hmm.LOOP] = hmm.xt[Hmm.XTN][Hmm.LOOP];
					command.newState.modified[0] = actionTarget.modified;
					commandHistory.add(command);
					layout.setHmmValue(hmm);
					layout.generateShape(layout.size);
					repaint();
				}
				return;
			} else if (actionTarget == layout.bNode) {
				tp = new Probabilities(hmm.M + 1);
				tp.add("D1", hmm.tbd1);
				for (int i = 1; i <= hmm.M; i++) {
					tp.add("M" + i, hmm.begin[i]);
				}
				JModifyNodeDialog dialog = new JModifyNodeDialog("B", tp, null);
				dialog.setVisible(true);
				if (dialog.action == JOptionPane.OK_OPTION) {
					HmmCommand command = new HmmCommand(hmm, layout, null);
					command.target = HmmCommand.Target.BNode;
					command.commandType = HmmCommand.CommandType.Modify;
					command.oldState.tbd1 = hmm.tbd1;
					command.oldState.beginArray = hmm.begin.clone();
					command.oldState.modified[0] = actionTarget.modified;

					actionTarget.modified = true;
					hmm.modified = true;
					hmm.tbd1 = tp.probs[0];
					for (int i = 1; i <= hmm.M; i++)
						hmm.begin[i] = tp.probs[i];

					command.newState.tbd1 = hmm.tbd1;
					command.newState.beginArray = hmm.begin.clone();
					command.newState.modified[0] = actionTarget.modified;
					commandHistory.add(command);

					layout.setHmmValue(hmm);
					layout.generateShape(layout.size);
					repaint();
				}
				return;
			} else if (actionTarget == layout.eNode) {
				tp = new Probabilities(2);
				tp.add("C", hmm.xt[Hmm.XTE][Hmm.MOVE]);
				tp.add("J", hmm.xt[Hmm.XTE][Hmm.LOOP]);
				JModifyNodeDialog dialog = new JModifyNodeDialog("E", tp, null);
				dialog.setVisible(true);
				if (dialog.action == JOptionPane.OK_OPTION) {
					HmmCommand command = new HmmCommand(hmm, layout, null);
					command.target = HmmCommand.Target.ENode;
					command.commandType = HmmCommand.CommandType.Modify;
					command.oldState.xt[Hmm.XTE][Hmm.MOVE] = hmm.xt[Hmm.XTE][Hmm.MOVE];
					command.oldState.xt[Hmm.XTE][Hmm.LOOP] = hmm.xt[Hmm.XTE][Hmm.LOOP];
					command.oldState.modified[0] = actionTarget.modified;

					actionTarget.modified = true;
					hmm.modified = true;
					hmm.xt[Hmm.XTE][Hmm.MOVE] = tp.probs[0];
					hmm.xt[Hmm.XTE][Hmm.LOOP] = tp.probs[1];

					command.newState.xt[Hmm.XTE][Hmm.MOVE] = hmm.xt[Hmm.XTE][Hmm.MOVE];
					command.newState.xt[Hmm.XTE][Hmm.LOOP] = hmm.xt[Hmm.XTE][Hmm.LOOP];
					command.newState.modified[0] = actionTarget.modified;
					commandHistory.add(command);

					layout.setHmmValue(hmm);
					layout.generateShape(layout.size);
					repaint();
				}
				return;
			} else if (actionTarget == layout.cNode) {
				tp = new Probabilities(2);
				tp.add("C", hmm.xt[Hmm.XTC][Hmm.MOVE]);
				tp.add("J", hmm.xt[Hmm.XTC][Hmm.LOOP]);
				JModifyNodeDialog dialog = new JModifyNodeDialog("C", tp, null);
				dialog.setVisible(true);
				if (dialog.action == JOptionPane.OK_OPTION) {
					HmmCommand command = new HmmCommand(hmm, layout, null);
					command.target = HmmCommand.Target.CNode;
					command.commandType = HmmCommand.CommandType.Modify;
					command.oldState.xt[Hmm.XTC][Hmm.MOVE] = hmm.xt[Hmm.XTC][Hmm.MOVE];
					command.oldState.xt[Hmm.XTC][Hmm.LOOP] = hmm.xt[Hmm.XTC][Hmm.LOOP];
					command.oldState.modified[0] = actionTarget.modified;

					actionTarget.modified = true;
					hmm.modified = true;
					hmm.xt[Hmm.XTC][Hmm.MOVE] = tp.probs[0];
					hmm.xt[Hmm.XTC][Hmm.LOOP] = tp.probs[1];

					command.newState.xt[Hmm.XTC][Hmm.MOVE] = hmm.xt[Hmm.XTC][Hmm.MOVE];
					command.newState.xt[Hmm.XTC][Hmm.LOOP] = hmm.xt[Hmm.XTC][Hmm.LOOP];
					command.newState.modified[0] = actionTarget.modified;
					commandHistory.add(command);

					layout.setHmmValue(hmm);
					layout.generateShape(layout.size);
					repaint();
				}
				return;
			} else if (actionTarget == layout.jNode) {
				tp = new Probabilities(2);
				tp.add("B", hmm.xt[Hmm.XTJ][Hmm.MOVE]);
				tp.add("J", hmm.xt[Hmm.XTJ][Hmm.LOOP]);
				JModifyNodeDialog dialog = new JModifyNodeDialog("J", tp, null);
				dialog.setVisible(true);
				if (dialog.action == JOptionPane.OK_OPTION) {
					HmmCommand command = new HmmCommand(hmm, layout, null);
					command.target = HmmCommand.Target.JNode;
					command.commandType = HmmCommand.CommandType.Modify;
					command.oldState.xt[Hmm.XTJ][Hmm.MOVE] = hmm.xt[Hmm.XTJ][Hmm.MOVE];
					command.oldState.xt[Hmm.XTJ][Hmm.LOOP] = hmm.xt[Hmm.XTJ][Hmm.LOOP];
					command.oldState.modified[0] = actionTarget.modified;

					actionTarget.modified = true;
					hmm.modified = true;
					hmm.xt[Hmm.XTJ][Hmm.MOVE] = tp.probs[0];
					hmm.xt[Hmm.XTJ][Hmm.LOOP] = tp.probs[1];

					command.newState.xt[Hmm.XTJ][Hmm.MOVE] = hmm.xt[Hmm.XTJ][Hmm.MOVE];
					command.newState.xt[Hmm.XTJ][Hmm.LOOP] = hmm.xt[Hmm.XTJ][Hmm.LOOP];
					command.newState.modified[0] = actionTarget.modified;
					commandHistory.add(command);

					layout.setHmmValue(hmm);
					layout.generateShape(layout.size);
					repaint();
				}
				return;
			}
		} else if (e.getActionCommand().equals("Delete")) {
			for (int i = 1; i <= hmm.M; i++) {
				if (actionTarget == layout.mNodes[i]) {
					int ret = JOptionPane.showConfirmDialog(null,
							"Do you want Remove State " + i + "?");
					if (ret == JOptionPane.YES_OPTION) {
						HmmCommand command = new HmmCommand(hmm, layout, null);
						command.target = HmmCommand.Target.ProfileHmm;
						command.commandType = HmmCommand.CommandType.Delete;
						command.position = i;
						hmm.deleteState(i, command);
						layout.deleteState(i, command);
						layout.setHmmValue(hmm);
						layout.generateShape(layout.size);
						clearSelection();
						repaint();
						commandHistory.add(command);
					}
				}
			}
		} else if (e.getActionCommand().equals("Before")) {
			for (int i = 1; i <= hmm.M; i++) {
				if (actionTarget == layout.mNodes[i]) {
					insertState(i, false);
					return;
				}
			}
		} else if (e.getActionCommand().equals("After")) {
			for (int i = 1; i <= hmm.M; i++) {
				if (actionTarget == layout.mNodes[i]) {
					insertState(i + 1, true);
					return;
				}
			}
		}
	}

	/**
	 * Insert state
	 * @param i position of the insertion
	 * @param before whether copy data from previous or next node
	 */
	void insertState(int i, boolean before) {
		HmmCommand command = new HmmCommand(hmm, layout, null);
		command.target = HmmCommand.Target.ProfileHmm;
		command.commandType = HmmCommand.CommandType.Insert;
		command.position = i;
		command.before = before;

		hmm.insertState(i, before, null);
		layout.insertState(i, before, null);
		layout.setHmmValue(hmm);
		layout.generateShape(layout.size);
		clearSelection();
		repaint();
		commandHistory.add(command);
	}
}
