#include "CurveControl.h"

#include "Geometry.h"
#include "Log.h"
#include "Panel.h"
#include "ShaderProgram.h"
#include "Window.h"

#include <imgui.h>
#include <memory>
#include <glm/glm.hpp>

class CurveEditorCallBack : public CallbackInterface {
public:
	CurveEditorCallBack() {}

	virtual void keyCallback(int key, int scancode, int action, int mods) override {
		Log::info("KeyCallback: key={}, action={}", key, action);
	}

	virtual void mouseButtonCallback(int button, int action, int mods) override {
		Log::info("MouseButtonCallback: button={}, action={}", button, action);

		if (button == GLFW_MOUSE_BUTTON_1 && action == GLFW_PRESS) {
			wasClicked = true;
		}
	}

	virtual void cursorPosCallback(double xpos, double ypos) override {
		Log::info("CursorPosCallback: xpos={}, ypos={}", xpos, ypos);
	}

	virtual void scrollCallback(double xoffset, double yoffset) override {
		Log::info("ScrollCallback: xoffset={}, yoffset={}", xoffset, yoffset);
	}

	virtual void windowSizeCallback(int width, int height) override {
		Log::info("WindowSizeCallback: width={}, height={}", width, height);

		// **Important**, calls glViewport(0, 0, width, height);
		CallbackInterface::windowSizeCallback(width, height);
	}

	bool wasClicked = false;
};

// Can swap the callback instead of maintaining a state machine

class TurnTable3DViewerCallBack : public CallbackInterface {
public:
	TurnTable3DViewerCallBack() {}

	virtual void keyCallback(int key, int scancode, int action, int mods) {}
	virtual void mouseButtonCallback(int button, int action, int mods) {}
	virtual void cursorPosCallback(double xpos, double ypos) {}
	virtual void scrollCallback(double xoffset, double yoffset) {}
	virtual void windowSizeCallback(int width, int height) {

		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width, height);
	}
};

class CurveEditorPanelRenderer : public PanelRendererInterface {
public:
	CurveEditorPanelRenderer()
		: inputText(""), buttonClickCount(0), pointSize(5.0f), dragValue(0.0f),
		inputValue(0.0f), checkboxValue(false), comboSelection(0) {
		// Initialize options for the combo box
		options[0] = "Option 1";
		options[1] = "Option 2";
		options[2] = "Option 3";

		// Initialize color (white by default)
		colorValue[0] = 1.0f; // R
		colorValue[1] = 1.0f; // G
		colorValue[2] = 1.0f; // B
	}

	virtual void render() override {
		// Color selector
		ImGui::ColorEdit3("Select Background Color",
			colorValue); // RGB color selector
		ImGui::Text("Selected Color: R: %.3f, G: %.3f, B: %.3f", colorValue[0],
			colorValue[1], colorValue[2]);

		// Text input
		ImGui::InputText("Input Text", inputText, IM_ARRAYSIZE(inputText));

		// Display the input text
		ImGui::Text("You entered: %s", inputText);

		// Button
		if (ImGui::Button("Click Me")) {
			buttonClickCount++;
		}
		ImGui::Text("Button clicked %d times", buttonClickCount);

		// Scrollable block
		ImGui::TextWrapped("Scrollable Block:");
		ImGui::BeginChild("ScrollableChild", ImVec2(0, 100),
			true); // Create a scrollable child
		for (int i = 0; i < 20; i++) {
			ImGui::Text("Item %d", i);
		}
		ImGui::EndChild();

		// Float slider
		ImGui::SliderFloat("Float Slider", &pointSize, 5.0f, 100.0f,
			"Point Size: %.3f");

		// Float drag
		ImGui::DragFloat("Float Drag", &dragValue, 0.1f, 0.0f, 100.0f,
			"Drag Value: %.3f");

		// Float input
		ImGui::InputFloat("Float Input", &inputValue, 0.1f, 1.0f,
			"Input Value: %.3f");

		// Checkbox
		ImGui::Checkbox("Enable Feature", &checkboxValue);
		ImGui::Text("Feature Enabled: %s", checkboxValue ? "Yes" : "No");

		// Combo box
		ImGui::Combo("Select an Option", &comboSelection, options,
			IM_ARRAYSIZE(options));
		ImGui::Text("Selected: %s", options[comboSelection]);

		// Displaying current values
		ImGui::Text("Point Size: %.3f", pointSize);
		ImGui::Text("Drag Value: %.3f", dragValue);
		ImGui::Text("Input Value: %.3f", inputValue);
	}

	glm::vec3 getColor() const {
		return glm::vec3(colorValue[0], colorValue[1], colorValue[2]);
	}

	float getPointSize() const { return pointSize; }

private:
	float colorValue[3];    // Array for RGB color values
	char inputText[256];    // Buffer for input text
	int buttonClickCount;   // Count button clicks
	float pointSize;        // Value for float slider
	float dragValue;        // Value for drag input
	float inputValue;       // Value for float input
	bool checkboxValue;     // Value for checkbox
	int comboSelection;     // Index of selected option in combo box
	const char* options[3]; // Options for the combo box
};

CurveControl::CurveControl(Window& window)
	: mShader("shaders/test.vert", "shaders/test.frag"),
	mPanel(window.getGLFWwindow()) {
	mCurveControls = std::make_shared<CurveEditorCallBack>();
	m3DCameraControls = std::make_shared<TurnTable3DViewerCallBack>();

	mPanelRenderer = std::make_shared<CurveEditorPanelRenderer>();

	// Set callback to window
	window.setCallbacks(mCurveControls);
	// Can swap the callback instead of maintaining a state machine
	// window.setCallbacks(turn_table_3D_viewer_callback);

	mPanel.setPanelRenderer(mPanelRenderer);

	mCurveGeometry = GenerateInitialGeometry();
	mGPUGeometry.setVerts(mCurveGeometry.verts);
	mGPUGeometry.setCols(mCurveGeometry.cols);

	// Using two different buffers for the control points and the lines themselves
	// makes it easier to highlight the selected point
	mPointGPUGeometry.setVerts(mCurveGeometry.verts);
	mPointGPUGeometry.setCols(
		std::vector<glm::vec3>(mCurveGeometry.verts.size(), { 1.f, 0.f, 0.f }));

	// Setup OpenGL state that stays constant for most of the program
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_FRAMEBUFFER_SRGB);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void CurveControl::DrawGeometry() {
	glm::vec3 backgroundColor = mPanelRenderer->getColor();

	glClearColor(backgroundColor.r, backgroundColor.g, backgroundColor.b, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	mShader.use();

	// Render the control points
	glPointSize(mPanelRenderer->getPointSize());
	mGPUGeometry.bind();
	glDrawArrays(GL_POINTS, 0, mCurveGeometry.verts.size());

	// Render the curve that connects the control points
	mPointGPUGeometry.bind();
	glDrawArrays(GL_LINE_STRIP, 0, mCurveGeometry.verts.size());

	// disable sRGB for things like imgui
	glDisable(GL_FRAMEBUFFER_SRGB);
	mPanel.render();
}

void CurveControl::Update() {
	// Use this function to process logic and update things based on user inputs
	//Example: generate a new control point

	if (mCurveControls->wasClicked) {
		Log::debug("Insert or select a control point based on the position clicked");
		mCurveControls->wasClicked = false;
	}
}

// Generate some initial points to show what the rendering system is doing at
// the moment
CPU_Geometry CurveControl::GenerateInitialGeometry() {
	std::vector<glm::vec3> cp_positions_vector = {
	  {-.5f, -.5f, 0.f}, {.5f, -.5f, 0.f}, {.5f, .5f, 0.f}, {-.5f, .5f, 0.f} };
	glm::vec3 cp_line_colour = { 0.f, 1.f, 0.f };

	CPU_Geometry cp_point_cpu;
	cp_point_cpu.verts = cp_positions_vector;
	cp_point_cpu.cols =
		std::vector<glm::vec3>(cp_point_cpu.verts.size(), cp_line_colour);

	return cp_point_cpu;
}
