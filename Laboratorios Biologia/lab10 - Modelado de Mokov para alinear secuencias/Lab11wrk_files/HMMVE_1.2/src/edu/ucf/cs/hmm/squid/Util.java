package edu.ucf.cs.hmm.squid;

public class Util {
	public static boolean isgap(char c) {
		if ((c) == ' ' || (c) == '.' || (c) == '_' || (c) == '-' || (c) == '~')
			return true;
		return false;
	}

	public static boolean isalpha(char c) {
		if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'))
			return true;
		return false;
	}
	public static boolean isdigit(char c) {
		if (c >= '0' && c <= '9')
			return true;
		return false;
	}
	public static boolean isspace(char c) {
		if (c == ' ' || c == '\t' || c == '\n')
			return true;
		return false;
	}
	public static boolean islower(char c) {
		if (c>='a'&&c<='z')
			return true;
		return false;
	}
	public static boolean isalnum(char c) {
		if ((c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9'))
			return true;
		return false;
	}
}
