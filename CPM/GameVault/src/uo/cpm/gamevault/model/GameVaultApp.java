package uo.cpm.gamevault.model;

public class GameVaultApp {
	
	private GameVault vault;
	
	public void start() {
		//TODO manejo de excepciones de entrada/salida 
		vault=new GameVault();
		vault.loadGenres();
		vault.loadGames();
		vault.loadLibrary();
		vault.run();
	}
}
