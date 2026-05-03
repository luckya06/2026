package uo.cpm.gamevault.model;

import java.util.List;
import uo.cpm.util.gamevault.file.*;

public class GameVault {
	List<Genre> genres;
	List<Game> games;

	public void loadGenres() {
		FileUtil.getLinesFromFile("files/generos.dat");
		
	}

	public void loadGames() {
		FileUtil.getLinesFromFile("files/juegos.dat");
		
	}

	public void loadLibrary() {
		// TODO Auto-generated method stub
		
	}

	public void run() {
		// TODO Auto-generated method stub
		
	}
	
	

}
