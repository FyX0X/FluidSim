#include <LittleEngine/little_engine.h>
#include "game.h"


int main()
{

	game::Game gameInstance;


	LittleEngine::EngineConfig config;
	config.title = "Template Test Game";
	config.iconPath = RESOURCES_PATH "awesomeface.png"; // Set the path to your icon if you have one

	LittleEngine::Initialize(config);	
	LittleEngine::GetWindow()->SetWindowResizeCallback([&](int w, int h) { gameInstance.OnWindowSizeChange(w, h); });

	gameInstance.Initialize();

	LittleEngine::Run(
		[&](float dt) { gameInstance.Update(dt); },
		[&]() { gameInstance.Render(); }
	);

	LittleEngine::Shutdown();

	return 0;
}