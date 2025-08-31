#pragma once
#include <LittleEngine/little_engine.h>

#include "gameData.h"


namespace game
{

	class Game
	{

	public:
		Game() {};
		~Game() { Shutdown(); };

		Game(Game& other) = delete;
		Game(Game&& other) = delete;
		Game operator=(Game other) = delete;
		Game operator=(Game& other) = delete;
		Game operator=(Game&& other) = delete;

		bool Initialize();
		void Shutdown();

		void Update(float dt);
		void Render();

		// must always be set correctly
		void OnWindowSizeChange(int w, int h);


	private:

		void InitializeEngine();
		void InitializeResources();
		void InitializeInput();
		void InitializeUI();
		void InitializeScene();
		void InitializeLight();



		void ResizeFBOs();

		void BlurLightTexture(LittleEngine::Graphics::RenderTarget& lightFBO, int passes, LittleEngine::Graphics::Shader& shader);
		


		



		GameData m_data;
		std::unique_ptr<LittleEngine::Graphics::Renderer> m_renderer;
		std::unique_ptr<LittleEngine::Audio::AudioSystem> m_audioSystem;
		std::unique_ptr<LittleEngine::UI::UISystem> m_uiSystem; // UI system for handling UI elements and contexts
		std::unique_ptr<LittleEngine::Graphics::LightSystem> m_lightSystem; // light system for rendering lights and shadows

		// temporary

		bool outlineMode = false;
		


		float speed = 10.f;


		bool w = false;
		bool v = true;
		// fps:
		float delta = 0;
		std::vector<float> fpsHistory;
		std::vector<float> distHistory;
		const int historySize = 1000;


		LittleEngine::Graphics::Camera sceneCamera = {};


	};





}