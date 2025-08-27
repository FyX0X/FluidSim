#include "game.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#if ENABLE_IMGUI == 1
#include "imgui.h"
#endif

#include <LittleEngine/little_engine.h>

#include <LittleEngine/Utils/logger.h>

// Temporary includes for glad and GLFW for keys
#include <glad/glad.h>

#ifdef USE_GLFW
#include <GLFW/glfw3.h>
#endif // USE_GLFW

#include <iostream>
#include <iomanip>
#include <sstream>



namespace game
{

	const unsigned int Game::world[] = { 2, 2, 2, 3, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

#pragma region Game Initialization/Unintialization

	bool Game::Initialize()
	{
		InitializeEngine();

		InitializeResources();

		InitializeInput();

		InitializeUI();

		InitializeScene();

		InitializeLight();


		




		// TESTING
		LittleEngine::Utils::Logger::Info("info");
		LittleEngine::Utils::Logger::Warning("warning");
		LittleEngine::Utils::Logger::Error("error");
		//LittleEngine::Utils::Logger::Critical("critical");

		glm::vec2 a = { 0, 0 };
		glm::vec2 b = { 1, 0 };
		glm::vec2 c = { 2, 0 };

		LittleEngine::Math::Edge e1 = { a, c };

		LittleEngine::Math::Edge e2 = { b, c };

		// triangle in CCW order

		std::cout << "Triangle signed area: " << LittleEngine::Math::TriangleSignedArea(a, b, c) << "\n";
		std::cout << "Three point orientation: " << LittleEngine::Math::ThreePointOrientation(a, b, c) << "\n";

		std::cout << "b on AC: " << LittleEngine::Math::PointOnSegment(b, e1) << "\n";
		std::cout << "a on BC: " << LittleEngine::Math::PointOnSegment(a, e2) << "\n";
		std::cout << "c on AC: " << LittleEngine::Math::PointOnSegment(c, e1) << "\n";
		std::cout << "(1, 0.1) on AC: " << LittleEngine::Math::PointOnSegment({ 1, 0.1f }, e1) << "\n";


		//loading the saved data. Loading an entire structure like this makes savind game data very easy.
		//platform::readEntireFile(RESOURCES_PATH "gameData.data", &gameData, sizeof(GameData));

		return true;
	}


	void Game::InitializeEngine()
	{
		m_renderer = std::make_unique<LittleEngine::Graphics::Renderer>();
		m_renderer->Initialize(sceneCamera, LittleEngine::GetWindowSize());

		m_lightSystem = std::make_unique<LittleEngine::Graphics::LightSystem>();
		m_lightSystem->Initialize(1000); // initialize light system with a maximum of 1000 shadow quads

		m_audioSystem = std::make_unique<LittleEngine::Audio::AudioSystem>();
		m_audioSystem->Initialize();

		m_uiSystem = std::make_unique<LittleEngine::UI::UISystem>();
		m_uiSystem->Initialize(LittleEngine::GetWindowSize()); // initialize UI system with the current window size
	}

	void Game::InitializeResources()
	{
		m_audioSystem->LoadSound(RESOURCES_PATH "test.wav", sound);
		sound.SetPitch(pitch);
		sound.SetVolume(volume);
		sound.SetSpatialization(spatialized);
		if (spatialized)
		{
			sound.SetMinDistance(minD);
			sound.SetMaxDistance(maxD);
			sound.SetMaxGain(maxG);
			sound.SetMinGain(minG);
			sound.SetRolloff(rolloff);
		}

		texture1.LoadFromFile(RESOURCES_PATH "test.jpg");
		texture2.LoadFromFile(RESOURCES_PATH "awesomeface.png", true, false, true);
		torch.LoadFromFile(RESOURCES_PATH "torch.png");
		minecraft_blocks.LoadFromFile(RESOURCES_PATH "minecraft_atlas.png");
		minecraft_atlas = LittleEngine::Graphics::TextureAtlas(minecraft_blocks, 16, 16);

		font.LoadFromTTF(RESOURCES_PATH "arial.ttf", 64.f);
	}

	void Game::InitializeInput()
	{

		// Register two separate axes
		LittleEngine::Input::RegisterAxis("vertical");
		LittleEngine::Input::RegisterAxis("horizontal");

		// Primary keys for each axis
		LittleEngine::Input::BindKeysToAxis(LittleEngine::Input::KeyCode::W, LittleEngine::Input::KeyCode::S, "vertical");
		LittleEngine::Input::BindKeysToAxis(LittleEngine::Input::KeyCode::D, LittleEngine::Input::KeyCode::A, "horizontal");

		// Duplicate (secondary) keys for each axis
		LittleEngine::Input::RegisterAxis("vertical_alt");
		LittleEngine::Input::RegisterAxis("horizontal_alt");

		LittleEngine::Input::BindKeysToAxis(LittleEngine::Input::KeyCode::Up, LittleEngine::Input::KeyCode::Down, "vertical_alt");
		LittleEngine::Input::BindKeysToAxis(LittleEngine::Input::KeyCode::Right, LittleEngine::Input::KeyCode::Left, "horizontal_alt");

		// Duplicate (third) keys for each axis
		LittleEngine::Input::RegisterAxis("vertical3");
		LittleEngine::Input::RegisterAxis("horizontal3");

		LittleEngine::Input::BindKeysToAxis(LittleEngine::Input::KeyCode::I, LittleEngine::Input::KeyCode::K, "vertical3");
		LittleEngine::Input::BindKeysToAxis(LittleEngine::Input::KeyCode::L, LittleEngine::Input::KeyCode::J, "horizontal3");


#pragma region Command definitions

		class ZoomCommand : public LittleEngine::Input::Command {
			float& zoom;
		public:
			ZoomCommand(float& zoom) : zoom(zoom) {}
			std::string GetName() const override { return "Zoom"; }

			void OnPress() override { zoom *= 2.f; }
			void OnHold() override {
				zoom *= 0.99f;
				if (zoom <= 1)
					zoom = 1.f;
			}
			void OnRelease() override {}
		};

		class JumpCommand : public LittleEngine::Input::Command {
			glm::vec2& pos;
		public:
			JumpCommand(glm::vec2& pos) : pos(pos) {}
			std::string GetName() const override { return "Jump"; }

			void OnPress() override { pos.y += 2.f; }
			void OnRelease() override { pos.y -= 2.f; }
			void OnHold() override {}
		};

		class ScreenshotCommand : public LittleEngine::Input::Command {
			LittleEngine::Graphics::Renderer* renderer;
		public:
			ScreenshotCommand(LittleEngine::Graphics::Renderer* r) : renderer(r) {}
			std::string GetName() const override { return "Screenshot"; }

			void OnPress() override { renderer->SaveScreenshot(); }
		};

		class SaveRenderTargetCommand : public LittleEngine::Input::Command {
			LittleEngine::Graphics::Renderer* renderer;
			LittleEngine::Graphics::RenderTarget& target;
		public:
			SaveRenderTargetCommand(LittleEngine::Graphics::Renderer* r,
				LittleEngine::Graphics::RenderTarget& t) : renderer(r), target(t) {
			}

			std::string GetName() const override { return "SaveRenderTarget"; }

			void OnPress() override {
				std::cout << "hello before: ";
				renderer->SaveScreenshot(&target);
				std::cout << "hello after\n ";
			}
		};


		class ColorCommand : public LittleEngine::Input::Command {
			LittleEngine::Graphics::Color& color;
		public:
			ColorCommand(LittleEngine::Graphics::Color& c) : color(c) {}
			std::string GetName() const override { return "Color"; }

			void OnPress() override {
				float f = (rand() % 1000) / 1000.f;
				color = f * LittleEngine::Graphics::Colors::Green;
				color.w = 1.f;
			}
			void OnRelease() override {}
			void OnHold() override {}
		};

		class SoundCommand : public LittleEngine::Input::Command {
			LittleEngine::Audio::Sound& sound;
		public:
			SoundCommand(LittleEngine::Audio::Sound& s) : sound(s) {}
			std::string GetName() const override { return "Sound"; }

			void OnPress() override {

				// play sound
				sound.Play();
			}
			void OnRelease() override {}
			void OnHold() override {}
		};

		class UpdatePointPosCommand : public LittleEngine::Input::Command {
			glm::vec2& pos;
			const LittleEngine::Graphics::Camera& camera;
		public:
			UpdatePointPosCommand(glm::vec2& p, LittleEngine::Graphics::Camera& c) : pos(p), camera(c) {}
			std::string GetName() const override { return "UpdatePointPos"; }

			void OnPress() override {
				// update position to mouse position
				glm::ivec2 screenPos = LittleEngine::Input::GetMousePosition();
				glm::vec2 ndc;
				ndc.x = (2.0f * screenPos.x) / camera.viewportSize.x - 1.0f;
				ndc.y = ((2.0f * screenPos.y) / camera.viewportSize.y - 1.0f);

				// Homogeneous clip space (z=0, w=1 for 2D)
				glm::vec4 clipPos(ndc, 0.0f, 1.0f);

				// Compute inverse view-projection
				glm::mat4 invViewProj = glm::inverse(camera.GetProjectionMatrix() * camera.GetViewMatrix());

				// Transform to world space
				glm::vec4 worldPos = invViewProj * clipPos;

				pos.x = worldPos.x / worldPos.w; // divide by w to get correct position
				pos.y = worldPos.y / worldPos.w; // divide by w to get correct position

			}
		};

		class AddPointCommand : public LittleEngine::Input::Command {
			LittleEngine::Math::Polygon& poly;
			const LittleEngine::Graphics::Camera& camera;
		public:
			AddPointCommand(LittleEngine::Math::Polygon& p, LittleEngine::Graphics::Camera& c) : poly(p), camera(c) {}
			std::string GetName() const override { return "AddPoint"; }
			void OnPress() override {
				// add point to polygon at mouse position
				glm::ivec2 screenPos = LittleEngine::Input::GetMousePosition();
				glm::vec2 ndc;
				ndc.x = (2.0f * screenPos.x) / camera.viewportSize.x - 1.0f;
				ndc.y = ((2.0f * screenPos.y) / camera.viewportSize.y - 1.0f);
				// Homogeneous clip space (z=0, w=1 for 2D)
				glm::vec4 clipPos(ndc, 0.0f, 1.0f);
				// Compute inverse view-projection
				glm::mat4 invViewProj = glm::inverse(camera.GetProjectionMatrix() * camera.GetViewMatrix());
				// Transform to world space
				glm::vec4 worldPos = invViewProj * clipPos;
				poly.vertices.push_back({ worldPos.x / worldPos.w, worldPos.y / worldPos.w });
			}
		};

		class ResetPolygonCommand : public LittleEngine::Input::Command {
			LittleEngine::Math::Polygon& poly;
		public:
			ResetPolygonCommand(LittleEngine::Math::Polygon& p) : poly(p) {}
			std::string GetName() const override { return "ResetPolygon"; }
			void OnPress() override {
				poly.vertices.clear();
			}
		};

#pragma endregion


		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::U, std::make_unique<UpdatePointPosCommand>(m_data.A, sceneCamera));
		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::I, std::make_unique<UpdatePointPosCommand>(m_data.B, sceneCamera));
		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::O, std::make_unique<UpdatePointPosCommand>(m_data.C, sceneCamera));
		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::P, std::make_unique<UpdatePointPosCommand>(m_data.D, sceneCamera));

		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::T, std::make_unique<AddPointCommand>(polygon, sceneCamera));
		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::R, std::make_unique<ResetPolygonCommand>(polygon));

		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::Space, std::make_unique<ColorCommand>(m_data.color));
		//LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::SPACE, std::make_unique<JumpCommand>(m_data.rectPos));
		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::Space, std::make_unique<SoundCommand>(sound));
		LittleEngine::Input::BindMouseButtonToCommand(LittleEngine::Input::MouseButton::Left, std::make_unique<ColorCommand>(m_data.color));
		//LittleEngine::Input::BindMouseButtonToCommand(LittleEngine::Input::MouseButton::Left, std::make_unique<ZoomCommand>(m_data.zoom));
		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::F11, std::make_unique<ScreenshotCommand>(m_renderer.get()));
		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::F10, std::make_unique<SaveRenderTargetCommand>(m_renderer.get(), target));

	}

	void Game::InitializeUI()
	{	
		LittleEngine::UI::UIContext* hud = m_uiSystem->CreateContext("HUD");
		LittleEngine::UI::UIContext* menu = m_uiSystem->CreateContext("Menu");


		// add elements
		glm::ivec2 size = LittleEngine::GetWindowSize();

		LittleEngine::UI::UIElement* label;
		label = hud->AddElement(std::make_unique<LittleEngine::UI::UILabel>(glm::vec2{5, size.y * 0.95}, "Hello World!", 32.f));
		
		std::unique_ptr<LittleEngine::UI::UIButton> button = std::make_unique<LittleEngine::UI::UIButton>(
			glm::vec2{ size.x / 2 - 50, size.x * 0.6f }, 
			glm::vec2{100, 50}, "Button");
		
		button->SetOnClickCallback([&]() {
			std::cout << "Button clicked!\n";
			m_uiSystem->ToggleContext("Menu");	// toggle HUD context
		});
		hud->AddElement(std::move(button));
		
		menu->AddElement(std::make_unique<LittleEngine::UI::UILabel>(glm::vec2{ size.x / 2 - 30, size.y * 0.9f }, "Menu", 32.f));

		std::unique_ptr<LittleEngine::UI::UICheckbox> cb = std::make_unique<LittleEngine::UI::UICheckbox>(
			glm::vec2{ size.x / 2 - 50, size.x * 0.4f }, 
			40.f, "Checkbox", false);

		cb->SetOnToggleCallback([&](bool state) {
			std::cout << "Checkbox state: " << (state ? "Checked" : "Unchecked") << "\n";
			});
		
		menu->AddElement(std::move(cb));

		m_uiSystem->ShowContext("HUD");	// show HUD context by default
	}

	void Game::InitializeScene()
	{
		ResizeFBOs();

		// Creating tilemap
		//				0		1		2		3
		tileIDs = { {3, 15}, {2, 15}, {1, 15}, {0, 15} };

		tilemap.SetTileSetTexture(minecraft_blocks, minecraft_atlas);
		//tilemap.SetTileSize(1.f);
		tilemap.SetTileSetAtlasKey(tileIDs);
		tilemap.SetMap(world, 10, 10, { 0, 0 });

		// create custom texture with render target

		sceneCamera.centered = true;	// center camera on screen
		sceneCamera.viewportSize = LittleEngine::GetWindowSize();

		target.Create(100, 100);
		m_renderer->SetRenderTarget(&target);

		LittleEngine::Graphics::Camera tempCamera = sceneCamera;
		tempCamera.centered = false;	// do not center camera on screen

		tempCamera.viewportSize = (glm::vec2)target.GetSize();
		m_renderer->SetCamera(tempCamera);

		m_renderer->DrawRect({ 0, 0, 100, 100 }, texture2);
		m_renderer->DrawString("target", { 10, 100 }, LittleEngine::Graphics::Colors::White, 10);

		m_renderer->Flush();
		m_renderer->SetRenderTarget();

		m_renderer->SetCamera(sceneCamera);



		for (int i = 0; i < 100; i++)
		{
			glm::vec4 r = { 0, 0, 0.8f, 0.8f };
			r.x = (i % 10) - 5;     // x goes from -5 to +4 across columns
			r.y = (i / 10) - 5;     // y goes from -5 to +4 across rows
			rect.push_back(r);
			if (i >= 90)
			{
				rect_uv.push_back(minecraft_atlas.GetUV(3, 15));
			}
			else
			{
				rect_uv.push_back(minecraft_atlas.GetUV(2, 15));
			}
		}
	}

	void Game::InitializeLight()
	{

		blurShader.Create(RESOURCES_PATH "fullscreen_quad.vert", RESOURCES_PATH "blur.frag", true);
		//lightSceneMergingShader.Create(RESOURCES_PATH "fullscreen_quad.vert", RESOURCES_PATH "merge_light_scene.frag", true);

		std::vector<glm::vec2> vertices = {
			{ -1.f, -1.f },
			{ 1.f, -1.f },
			{ 1.f, 1.f },
			{ -1.f, 1.f }
		};
		obstacles.push_back(m_lightSystem->CreateObstacle(vertices));

		vertices = {

			{ -3.f, -3.f },
			{ -2.f, -3.f },
			{ -1.f, -2.f }
		};
		obstacles.push_back(m_lightSystem->CreateObstacle(vertices));

		vertices = {
			{ -4.f, 3.f },
			{ -3.f, 4.f },
			{ -4.f, 5.f }
		};
		obstacles.push_back(m_lightSystem->CreateObstacle(vertices));

		lightSources.push_back(m_lightSystem->CreateLightSource({ 0.f, 0.f }, { 1.f, 1.f, 1.f }, 1.f, 10.f));
		lightSources.push_back(m_lightSystem->CreateLightSource({ 3.f, 3.f }, { 0.f, 1.f, 1.f }, 1.f, 4.f));
		lightSources.push_back(m_lightSystem->CreateLightSource({ -3.f, -3.f }, { 1.f, 0.f, .5f }, 2.f, 15.f));


	}

	void Game::Shutdown()
	{
		m_renderer->Shutdown();
		m_audioSystem->Shutdown();
		sound.Shutdown();

		//saved the data.
		// platform::writeEntireFile(RESOURCES_PATH "gameData.data", &gameData, sizeof(GameData));
	}


#pragma endregion

#pragma region Mainloop

	void Game::Update(float dt)
	{
		delta = dt;


		// check axis input.

		float vert = LittleEngine::Input::GetAxis("vertical");
		float hori = LittleEngine::Input::GetAxis("horizontal");

		glm::vec2 move = { hori, vert };
		float length = glm::length(move);
		if (length > 1)
			move /= length;

		m_data.rectPos += move * speed * dt;

		vert = LittleEngine::Input::GetAxis("vertical_alt");
		hori = LittleEngine::Input::GetAxis("horizontal_alt");

		move = { hori, vert };
		length = glm::length(move);
		if (length > 1)
			move /= length;

		m_data.pos2 += move * 10.f * dt;

		vert = LittleEngine::Input::GetAxis("vertical3");
		hori = LittleEngine::Input::GetAxis("horizontal3");

		move = { hori, vert };
		length = glm::length(move);
		if (length > 1)
			move /= length;

		lightSources[2]->position += move * 10.f * dt;



		//m_renderer->camera.Follow(m_data.rectPos, dt, cameraFollowSpeed, maxDist, minDist);
		sceneCamera.FollowSpring(m_data.rectPos, dt, maxDist, cameraFollowSpeed);
		sceneCamera.zoom = m_data.zoom;

		m_audioSystem->SetListenerPosition(m_data.rectPos.x, m_data.rectPos.y);


		m_uiSystem->Update();

	}


	void Game::Render()
	{
#pragma region Game Rendering




#pragma region Scene Render
		// render scene to fbo
		m_renderer->BeginFrame();
		m_renderer->SetRenderTarget();
		m_renderer->SetRenderTarget(&sceneFBO);
		m_renderer->Clear();


		// green block from (-10, -10) to (0 0)
		m_renderer->DrawRect({ -10, -10, 10 , 10 }, m_data.color);

		// font block from (0, 0) to (15 15)
		//m_renderer->DrawRect({ 0, 0, 15 , 15 }, font.GetTexture());


		//m_renderer->DrawRect({ -15, 10, 15 , 15 }, texture2);



		for (int i = 0; i < rect.size(); i++)
		{
			m_renderer->DrawRect(rect[i], minecraft_blocks, color, rect_uv[i]);
			//m_renderer->DrawRect(rect[i], textures[i], color);
		}

		for (size_t i = 0; i < length; i++)
		{
			tilemap.Draw(m_renderer.get());		// each call draws the whole timeMap, only for benchmark purposes
		}

		m_renderer->DrawString("LittleEngine Template", { -.97, 4.97 }, font, LittleEngine::Graphics::Colors::Gray, scale);
		m_renderer->DrawString("LittleEngine Template", { -1, 5 }, font, LittleEngine::Graphics::Colors::White, scale);

		m_renderer->DrawString("Hello Arial", { 0, 0 }, font, LittleEngine::Graphics::Colors::White, scale);
		m_renderer->DrawString("Hello Default font", { 0, -3 }, LittleEngine::Graphics::Colors::White, scale);




		m_renderer->DrawRect(glm::vec4(m_data.pos2, 3.f, 3.f), target.GetTexture());


		LittleEngine::Math::Edge e1 = { m_data.A, m_data.B };
		LittleEngine::Math::Edge e2 = { m_data.C, m_data.D };

		LittleEngine::Graphics::Color c;
		if (LittleEngine::Math::SegmentsIntersect(e1, e2))
		{
			c = LittleEngine::Graphics::Colors::Green;
		}
		else
		{
			c = LittleEngine::Graphics::Colors::Red;
		}
		
		m_renderer->DrawLine(e1, 0.1f, c);
		m_renderer->DrawLine(e2, 0.1f, c);


		// draw polygon
		if (polygon.IsValid())
		{
			if (outlineMode)
			{
				m_renderer->DrawPolygonOutline(polygon, 0.1f, LittleEngine::Graphics::Colors::White);
			}
			else
			{
				m_renderer->DrawPolygon(polygon, LittleEngine::Graphics::Colors::White);
			}
		}




		m_renderer->Flush();

#pragma endregion


#pragma region Light Render

		m_lightSystem->RenderLighting(m_renderer.get(), &lightFBO, enableShadows);

#pragma endregion

		m_renderer->SetRenderTarget();
		m_renderer->MergeLightScene(lightFBO.GetTexture(), sceneFBO.GetTexture());



		m_renderer->shader.Use();


		// render UI;

		m_uiSystem->Render(m_renderer.get());

#pragma endregion



#pragma region ImGui Render
#if ENABLE_IMGUI == 1
 		float fps = 1.f / delta;
		// Your ImGui widgets here
		ImGui::Begin("Debug");
		ImGui::Text("FPS: %.2f", LittleEngine::GetFPS());
		ImGui::Text("QuadCount: %d", m_renderer->GetQuadCount());
		ImGui::Text("camera pos: %.1f, %.1f", sceneCamera.position.x, sceneCamera.position.y);
		ImGui::SliderFloat("Camera Zoom", &m_data.zoom, 0.1f, 100.f);
		ImGui::SliderFloat("light intensity", &lightIntensity, 0.1f, 100.f);
		ImGui::SliderInt("Blur passes", &blurPasses, 1, 20);
		if (ImGui::SliderInt("Downscale Factor", &downscaleFactor, 1, 20))
		{
			// resize the render targets
			ResizeFBOs();
		}
		ImGui::Checkbox("Enable Shadows", &enableShadows);
		ImGui::Checkbox("Outline Mode", &outlineMode);
		//ImGui::SliderFloat("Camera x", &m_data.rectPos.x, -50.f, 50.f);
		//ImGui::SliderFloat("Camera y", &m_data.rectPos.y, -50.f, 50.f);
		//ImGui::SliderFloat("Red", &color.x, 0.f, 1.f);
		//ImGui::SliderFloat("Green", &color.y, 0.f, 1.f);
		//ImGui::SliderFloat("Blue", &color.z, 0.f, 1.f);
		//ImGui::SliderFloat("Alpha", &color.w, 0.f, 1.f);
		// 
		// SOUND
		// 
		//if (ImGui::SliderFloat("min dist", &minD, 0.f, 20.f))
		//{
		//	sound.SetMinDistance(minD);
		//}
		//if (ImGui::SliderFloat("max dist", &maxD, 0.f, 20.f))
		//{
		//	sound.SetMaxDistance(maxD);
		//}		
		//if (ImGui::SliderFloat("min gain", &minG, 0.f, 1.f))
		//{
		//	sound.SetMinGain(minG);
		//}
		//if (ImGui::SliderFloat("max gain", &maxG, 0.f, 1.f))
		//{
		//	sound.SetMaxGain(maxG);
		//}
		//if (ImGui::SliderFloat("rolloff", &rolloff, 0.f, 10.f))
		//{
		//	sound.SetRolloff(rolloff);
		//}		
		// if (ImGui::SliderFloat("pitch", &pitch, 0.f, 10.f))
		//{
		//	sound.SetPitch(pitch);
		//}
		//if (ImGui::SliderFloat("volume", &volume, 0.f, 10.f))
		//{
		//	sound.SetVolume(volume);
		//}
		ImGui::SliderFloat("scale", &scale, 0.1f, 5.f);
		ImGui::SliderFloat("speed", &speed, 0.f, 100.f);
		ImGui::SliderFloat("camera follow speed", &cameraFollowSpeed, 0.f, 30.f);
		ImGui::SliderFloat("max camera dist", &maxDist, 0.f, 10.f);
		ImGui::SliderFloat("min camera dist", &minDist, 0.f, 1.f);

		ImGui::SliderInt("Tilemap Count", &length, 0, 200);

		if (ImGui::Checkbox("wireframe", &w))
		{
			m_renderer->SetWireframe(w);
		}
		if (ImGui::Checkbox("vsinc", &v))
		{
			LittleEngine::SetVsync(v);
		}
		if (ImGui::Checkbox("spatialized", &spatialized))
		{
			sound.SetSpatialization(spatialized);
		}

		// fps graph
		if (fpsHistory.size() >= historySize)
		{
			fpsHistory.erase(fpsHistory.begin());
		}
		fpsHistory.push_back(LittleEngine::GetFPS());

		float maxfps = 0;
		for (float f : fpsHistory)
		{
			if (f > maxfps)
				maxfps = f;
		}

		// Plot the FPS graph
		ImGui::PlotLines("FPS Graph", fpsHistory.data(), fpsHistory.size(), 0,
			nullptr, 0.0f, maxfps * 1.5f, ImVec2(0, 80));

		// camera distance
		if (distHistory.size() >= historySize)
		{
			distHistory.erase(distHistory.begin());
		}
		distHistory.push_back(sceneCamera.position.x - m_data.rectPos.x);

		float maxx = 0;
		for (float f : distHistory)
		{
			if (f > maxx)
				maxx = f;
			if (-f > maxx)
				maxx = -f;
		}

		// Plot the FPS graph
		ImGui::PlotLines("Delta x Graph", distHistory.data(), distHistory.size(), 0,
			nullptr, -maxx, maxx * 1.2f, ImVec2(0, 80));


		ImGui::End();

#endif
#pragma endregion



	}


#pragma endregion

#pragma region Callbacks

	void Game::OnWindowSizeChange(int w, int h)
	{
		m_renderer->UpdateWindowSize(w, h);
		m_uiSystem->UpdateWindowSize(w, h);
		
		LittleEngine::Input::UpdateWindowSize(w, h);

		// update the camera
		sceneCamera.viewportSize = glm::ivec2(w, h);

		ResizeFBOs();
	}


#pragma endregion

	void Game::ResizeFBOs()
	{
		// resize the render targets
		sceneFBO.Cleanup();
		lightFBO.Cleanup();

		sceneFBO.Create(LittleEngine::GetWindowSize().x, LittleEngine::GetWindowSize().y, GL_RGB);
		lightFBO.Create(LittleEngine::GetWindowSize().x / downscaleFactor, LittleEngine::GetWindowSize().y / downscaleFactor, GL_RGB16F);

	}
	void Game::BlurLightTexture(LittleEngine::Graphics::RenderTarget& lightFBO, int passes, LittleEngine::Graphics::Shader& shader)
	{
		bool horizontal = true;
		bool firstIteration = true;

		// 1. Create two ping-pong FBOs
		LittleEngine::Graphics::RenderTarget pingpongFBO[2];
		glm::ivec2 size = lightFBO.GetSize();

		for (int i = 0; i < 2; ++i)
			pingpongFBO[i].Create(size.x, size.y); // should use RGB16F if possible

		// 2. Blur Passes
		shader.Use();

		shader.SetInt("image", 0);

		for (int i = 0; i < passes * 2; ++i) {
			pingpongFBO[horizontal].Bind();
			shader.SetBool("horizontal", horizontal);
			shader.SetFloat("texelSize", horizontal ? (1.f / size.x) : (1.f / size.y));

			if (firstIteration)
				lightFBO.GetTexture().Bind(0);
			else
				pingpongFBO[!horizontal].GetTexture().Bind(0);


			m_renderer->FlushFullscreenQuad();

			horizontal = !horizontal;
			if (firstIteration)
				firstIteration = false;
		}

		// 3. Final pass: write blurred result back into lightFBO
		lightFBO.Bind();  // bind original light FBO as target
		shader.SetBool("horizontal", false); // doesn't matter, we just draw the final result
		shader.SetFloat("texelSize", 1.f/ size.y); // doesn't matter, we just draw the final result
		pingpongFBO[!horizontal].GetTexture().Bind(0); // final blurred texture

		shader.SetInt("image", 0);

		m_renderer->FlushFullscreenQuad();

		// 4. Cleanup temporary FBOs
		for (int i = 0; i < 2; ++i)
			pingpongFBO[i].Cleanup();

		lightFBO.Unbind(); // unbind the light FBO to restore default framebuffer

		m_renderer->shader.Use(); // restore default shader
	}

}

