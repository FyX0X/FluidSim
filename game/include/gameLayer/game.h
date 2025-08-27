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
		LittleEngine::Math::Polygon polygon = {};


		std::vector<LittleEngine::Math::Polygon*> obstacles;
		std::vector<LittleEngine::Graphics::LightSource*> lightSources;

		LittleEngine::Audio::Sound sound;
		float pitch = 1.f;
		float volume = 1.f;
		bool spatialized = false;
		float minD = .5f;
		float maxD = 10.f;
		float maxG = 1.f;
		float minG = .1f;
		float rolloff = 1.f;

		LittleEngine::Graphics::Texture texture1, texture2;
		LittleEngine::Graphics::Texture torch;
		LittleEngine::Graphics::Texture minecraft_blocks;
		LittleEngine::Graphics::TextureAtlas minecraft_atlas;
		std::vector<LittleEngine::Graphics::Texture> textures;
		std::vector<glm::vec4> rect;
		std::vector<glm::vec4> rect_uv;
		LittleEngine::Graphics::Color color = LittleEngine::Graphics::Colors::White;
		LittleEngine::Graphics::Font font;

		float scale = 1.f;

		float delta = 0;

		float speed = 10.f;


		int blurPasses = 5;
		int downscaleFactor = 2;
		float lightIntensity = 1.f;

		float cameraFollowSpeed = 30.f;
		float maxDist = 5.f;
		float minDist = 0.f;

		static const unsigned int world[];

		LittleEngine::Graphics::TilemapRenderer tilemap;
		std::vector<LittleEngine::Graphics::AtlasCoord> tileIDs;

		int length = 1;
		bool w = false;
		bool v = true;
		// fps:
		std::vector<float> fpsHistory;
		std::vector<float> distHistory;
		const int historySize = 1000;


		LittleEngine::Graphics::RenderTarget target = {};
		LittleEngine::Graphics::RenderTarget sceneFBO = {};
		LittleEngine::Graphics::RenderTarget lightFBO = {};

		LittleEngine::Graphics::Camera sceneCamera = {};


		LittleEngine::Graphics::Shader blurShader = {};
		//LittleEngine::Graphics::Shader lightSceneMergingShader = {};

		//LittleEngine::Graphics::Shader lightShader = {};
		//LittleEngine::Graphics::Shader shadowShader = {};


		//LittleEngine::Graphics::Shader fullscreenShader = {};
		//LittleEngine::Graphics::Shader fullscreenImageBlitShader = {};

		bool enableShadows = true;


		LittleEngine::UI::UICheckbox* cb_ptr = nullptr;

	};





}