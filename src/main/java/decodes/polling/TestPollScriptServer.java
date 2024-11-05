/*
 * $Id$
 * 
 * This software was written by Cove Software, LLC ("COVE") under contract
 * to Alberta Environment and Sustainable Resource Development (Alberta ESRD).
 * No warranty is provided or implied other than specific contractual terms 
 * between COVE and Alberta ESRD.
 *
 * Copyright 2014 Alberta Environment and Sustainable Resource Development.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package decodes.polling;

import java.io.File;
import java.io.IOException;
import java.net.Socket;
import java.util.Date;

import ilex.net.BasicServer;
import ilex.net.BasicSvrThread;
import ilex.util.Logger;

public class TestPollScriptServer extends BasicServer
{
	String scriptFileName = null;

	public TestPollScriptServer(int port, String scriptFileName)
		throws IllegalArgumentException, IOException
	{
		super(port);
		this.scriptFileName = scriptFileName;
	}

	@Override
	protected BasicSvrThread newSvrThread(Socket sock) throws IOException
	{
		System.out.println("New Client");
		return new TestPollScriptServerThread(this, sock);
	}
	
	/**
	 * Usage: TestPollScriptServer listeningPort scriptname
	 * Open listening socket on specified port. When a client connects spawn a
	 * TestPollScriptServerThread to execute the specified script with the client.
	 * 
	 * @param args Two arguments expected: listeningPort scriptname
	 */
	public static void main(String args[])
		throws Exception
	{
		int port = Integer.parseInt(args[0]);
		Logger.instance().setMinLogPriority(Logger.E_DEBUG3);
		TestPollScriptServer tts = new TestPollScriptServer(port, args[1]);
		tts.listen();
	}
}

class TestPollScriptServerThread
	extends BasicSvrThread
{
	protected TestPollScriptServerThread(BasicServer parent, Socket socket)
	{
		super(parent, socket);
	}

	@Override
	protected void serviceClient()
	{
		PollScriptProtocol prot = new PollScriptProtocol();
		try
		{
			File scriptFile = new File(((TestPollScriptServer)parent).scriptFileName);
			System.out.println("New client detected. Reading script '" + scriptFile.getPath() + "'");
			prot.readScript(scriptFile);
			
			System.out.println("Constructing ioPort...");
			// mock up an IOPort to execute the script
			IOPort ioPort = new IOPort(null, 0, null);
			ioPort.setIn(getSocket().getInputStream());
			ioPort.setOut(getSocket().getOutputStream());
			
			System.out.println("Executing script...");
			prot.executeScript(ioPort, new Date(System.currentTimeMillis() - 3600000L));
		}
		catch (Exception ex)
		{
			System.err.println(ex);
			ex.printStackTrace(System.err);
		}
		disconnect();
	}
}