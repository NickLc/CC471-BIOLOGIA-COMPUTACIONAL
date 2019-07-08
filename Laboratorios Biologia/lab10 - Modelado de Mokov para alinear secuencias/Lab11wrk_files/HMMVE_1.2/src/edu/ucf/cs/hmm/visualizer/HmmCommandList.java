package edu.ucf.cs.hmm.visualizer;

import java.util.ArrayList;

import javax.swing.JMenuItem;

/**
 * List for all commands performed, to support undo and redo operation
 *
 */
public class HmmCommandList {
	JMenuItem undoMenuItem; // we need a reference to the menu item to
							// enable/disable it

	JMenuItem redoMenuItem;

	ArrayList<HmmCommand> commands;

	int currentPosition;

	HmmCommandList(JMenuItem _undoMenuItem, JMenuItem _redoMenuItem) {
		commands = new ArrayList<HmmCommand>();
		undoMenuItem = _undoMenuItem;
		redoMenuItem = _redoMenuItem;
		currentPosition = -1; // tail of the command list, currentPosition is
								// a valid command, it is the first command to
								// undo,
		// currentPosition+1 is the first command to redo
	}

	/**
	 * Undo one operation
	 * 
	 */
	HmmCommand.Target undo() {
		HmmCommand c = commands.get(currentPosition);
		c.undo();
		currentPosition--;
		if (currentPosition == -1)
			undoMenuItem.setEnabled(false);
		redoMenuItem.setEnabled(true);
		return c.target;
	}

	/**
	 * Redo one operation
	 */
	HmmCommand.Target redo() {
		HmmCommand c = commands.get(currentPosition + 1);
		c.redo();
		currentPosition++;
		if (currentPosition == commands.size() - 1)
			redoMenuItem.setEnabled(false);
		undoMenuItem.setEnabled(true);
		return c.target;
	}

	/**
	 * Add one operation to history list, in case of undo
	 * @param command
	 */
	void add(HmmCommand command) {
		while (currentPosition != commands.size() - 1) {
			commands.remove(commands.size() - 1);
		}
		commands.add(command);
		currentPosition++;
		redoMenuItem.setEnabled(false);
		undoMenuItem.setEnabled(true);
	}

	/**
	 * Clear operation history
	 *
	 */
	void clear() {
		currentPosition = -1;
		commands.clear();
		undoMenuItem.setEnabled(false);
		redoMenuItem.setEnabled(false);
	}
}
