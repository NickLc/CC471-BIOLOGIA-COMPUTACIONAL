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
 * Panel to show null model layout
 * 
 */
public class NullModelPanel extends JPanel implements MouseMotionListener,
		MouseListener, ActionListener {
	Hmm hmm;

	NullModelLayout layout;

	NodeState lingerShape = null;

	NodeState selectedState = null;

	Point selectedPosition = null;

	Point lastDragDrawPosition = null;

	StateChartLabel stateChartLabel;

	HmmCommandList commandHistory;

	JPopupMenu popupMenu;

	JScrollPane container;

	public NullModelPanel(NullModelLayout _layout, Hmm _hmm,
			StateChartLabel _stateChartLabel, HmmCommandList _commandHistory,
			JScrollPane _container) {
		layout = _layout;
		hmm = _hmm;
		commandHistory = _commandHistory;
		container = _container;
		stateChartLabel = _stateChartLabel;
		setPreferredSize(layout.getDimension());
		addMouseMotionListener(this);
		addMouseListener(this);
		setAutoscrolls(true);
		popupMenu = new JPopupMenu();
		JMenuItem menuItem = new JMenuItem("Modify Node...", KeyEvent.VK_E);
		menuItem.setActionCommand("Modify");
		menuItem.addActionListener(this);
		popupMenu.add(menuItem);
	}

	public void resize() {
		setPreferredSize(layout.getDimension());
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Graphics2D g2d = (Graphics2D) g;
		Rectangle2D.Float clip = new Rectangle2D.Float();
		clip.x = container.getViewport().getViewPosition().x;
		clip.y = container.getViewport().getViewPosition().y;
		clip.width = container.getViewport().getExtentSize().width;
		clip.height = container.getViewport().getExtentSize().height;
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
		if (hmm != null) {
			Graphics2D g2d = (Graphics2D) getGraphics();
			if (lingerShape != null && lingerShape.contains(e.getPoint()))
				return;
			if (lingerShape != null) {
				lingerShape.toggleSelected(g2d, false);
				lingerShape = null;
			}
			if (layout.gNode.contains(e.getPoint())) {
				lingerShape = layout.gNode;
				lingerShape.toggleSelected(g2d, true);
				stateChartLabel.setState("G", hmm.nullModel);
				stateChartLabel.repaint();
			}
		}
	}
	/**
	 * Double click or right click
	 */
	public void mouseClicked(MouseEvent e) {
		if (e.getButton() == MouseEvent.BUTTON1 && e.getClickCount() != 2
				|| e.getButton() == MouseEvent.BUTTON3
				&& e.getClickCount() != 1) {
			return;
		}
		ActionEvent ae = new ActionEvent(popupMenu,
				ActionEvent.ACTION_PERFORMED, "Modify");
		if (layout.gNode.contains(e.getPoint())) {
			if (e.getButton() == MouseEvent.BUTTON3)
				popupMenu.show(e.getComponent(), e.getX(), e.getY());
			else
				actionPerformed(ae);
			return;
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
		NodeState state = null;
		if (layout.gNode.contains(e.getPoint())) {
			state = layout.gNode;
		} else if (layout.fNode.contains(e.getPoint())) {
			state = layout.fNode;
		}
		if (state != null) {
			selectedState = state;
			selectedPosition = e.getPoint();
		} else {
			selectedState = null;
			selectedPosition = null;
		}
	}

	/**
	 * Release dragging
	 */
	public void mouseReleased(MouseEvent e) {
		if (lastDragDrawPosition != null) {
			layout.moveNode(selectedState, lastDragDrawPosition.x,
					lastDragDrawPosition.y);
			repaint();
		}
		selectedState = null;
		selectedPosition = null;
		lastDragDrawPosition = null;
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("Modify")) {
			Probabilities tp;
			Probabilities ep = new Probabilities(hmm.alphabets);
			tp = new Probabilities(2);
			tp.add("G", hmm.p1);
			tp.add("F", 1 - hmm.p1);
			ep.set(hmm.nullModel);
			JModifyNodeDialog dialog = new JModifyNodeDialog("G", tp, ep);
			dialog.setVisible(true);
			if (dialog.action == JOptionPane.OK_OPTION) {
				HmmCommand command;
				command = new HmmCommand(hmm, null, layout);
				command.target = HmmCommand.Target.NullHmm;
				command.commandType = HmmCommand.CommandType.Modify;
				command.oldState.nullModel = hmm.nullModel.clone();
				command.oldState.p1 = hmm.p1;
				command.oldState.modified[0] = layout.gNode.modified;

				hmm.nullModel = ep.probs;
				layout.gNode.modified = true;
				hmm.modified = true;
				hmm.p1 = tp.probs[0];
				layout.setHmmValue(hmm);

				command.newState.nullModel = hmm.nullModel.clone();
				command.newState.p1 = hmm.p1;
				command.newState.modified[0] = layout.gNode.modified;
				commandHistory.add(command);
				layout.generateShape(layout.size);
				repaint();
			}
		}
	}
}
