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

#pragma region Game Initialization/Unintialization

	bool Game::Initialize()
	{
		InitializeEngine();

		InitializeResources();

		InitializeInput();

		InitializeUI();

		InitializeScene();

		InitializeLight();



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



#pragma region Command definitions

		
		class ScreenshotCommand : public LittleEngine::Input::Command {
			LittleEngine::Graphics::Renderer* renderer;
		public:
			ScreenshotCommand(LittleEngine::Graphics::Renderer* r) : renderer(r) {}
			std::string GetName() const override { return "Screenshot"; }

			void OnPress() override { renderer->SaveScreenshot(); }
		};


		class AddSourceCommand : public LittleEngine::Input::Command {
			Fluid2D& fluid;
			const LittleEngine::Graphics::Camera& camera;
		public:
			AddSourceCommand(Fluid2D& f, LittleEngine::Graphics::Camera& c) : fluid(f), camera(c) {}
			std::string GetName() const override { return "AddSource"; }
			void OnHold() override {
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
				fluid.AddSource(static_cast<int>(worldPos.x / worldPos.w), static_cast<int>(worldPos.y / worldPos.w), 1.f);
			}
		};


#pragma endregion

		LittleEngine::Input::BindKeyToCommand(LittleEngine::Input::KeyCode::F11, std::make_unique<ScreenshotCommand>(m_renderer.get()));
		LittleEngine::Input::BindMouseButtonToCommand(LittleEngine::Input::MouseButton::Left, std::make_unique<AddSourceCommand>(m_data.fluid, sceneCamera));

	}

	void Game::InitializeUI()
	{

	}

	void Game::InitializeScene()
	{
		ResizeFBOs();

		m_renderer->SetRenderTarget(nullptr); // set to default framebuffer


		m_data.fluid.Initialize(20, 20, 1.f, 0.001f);


		sceneCamera.centered = true;	// center camera on screen
		sceneCamera.viewportSize = LittleEngine::GetWindowSize();

		m_renderer->SetCamera(sceneCamera);
		
	}

	void Game::InitializeLight()
	{

		
	}

	void Game::Shutdown()
	{
		m_renderer->Shutdown();
		m_audioSystem->Shutdown();

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

		m_data.cameraPos += move * speed * dt;




		//m_renderer->camera.Follow(m_data.rectPos, dt, cameraFollowSpeed, maxDist, minDist);
		sceneCamera.FollowSpring(m_data.cameraPos, dt, 5, 30);

		m_audioSystem->SetListenerPosition(m_data.cameraPos.x, m_data.cameraPos.y);


		m_uiSystem->Update();



		//m_data.fluid.UpdateTimeStep(dt);





	}


	void Game::Render()
	{
#pragma region Game Rendering




#pragma region Scene Render
		// render scene to fbo
		m_renderer->BeginFrame();
		m_renderer->Clear(LittleEngine::Graphics::Colors::Black);


		//m_data.fluid.Render(m_renderer.get());


		// test
		m_renderer->DrawRect(glm::vec4{ 0.f, 0.f, 1.f, 1.f });
		m_renderer->DrawRect({ -10, -10, 10 , 10 }, LittleEngine::Graphics::Colors::Magenta);


		std::cout << m_renderer->GetQuadCount() << " quads drawn.\n";
		std::cout << ((m_renderer->GetRenderTarget() == nullptr) ? "Default framebuffer\n" : "Custom framebuffer\n") << std::endl;
		std::cout << "size: " << LittleEngine::GetWindowSize().x << "x" << LittleEngine::GetWindowSize().y << std::endl;
		std::cout << "camera: " << sceneCamera.viewportSize.x << "x" << sceneCamera.viewportSize.y << std::endl;
		m_renderer->Flush();




		// render UI;

		//m_uiSystem->Render(m_renderer.get());

#pragma endregion



#pragma region ImGui Render
#if ENABLE_IMGUI == 1
 		float fps = 1.f / delta;
		// Your ImGui widgets here
		ImGui::Begin("Debug");
		ImGui::Text("FPS: %.2f", LittleEngine::GetFPS());
		ImGui::Text("QuadCount: %d", m_renderer->GetQuadCount());
		ImGui::Text("camera pos: %.1f, %.1f", sceneCamera.position.x, sceneCamera.position.y);
		ImGui::SliderFloat("Camera Zoom", &sceneCamera.zoom, 0.1f, 100.f);
		ImGui::Checkbox("Outline Mode", &outlineMode);
		ImGui::SliderFloat("speed", &speed, 0.f, 100.f);

		if (ImGui::Checkbox("wireframe", &w))
		{
			m_renderer->SetWireframe(w);
		}
		if (ImGui::Checkbox("vsinc", &v))
		{
			LittleEngine::SetVsync(v);
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
		distHistory.push_back(sceneCamera.position.x - m_data.cameraPos.x);

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

