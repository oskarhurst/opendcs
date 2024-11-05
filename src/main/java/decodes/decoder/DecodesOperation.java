/*
*	$Id$
*
*	$State$
*
*	$Log$
*	Revision 1.1  2008/04/04 18:21:01  cvs
*	Added legacy code to repository
*	
*	Revision 1.4  2004/08/31 16:31:19  mjmaloney
*	javadoc
*	
*	Revision 1.3  2001/08/19 19:33:21  mike
*	dev
*	
*	Revision 1.2  2001/05/06 22:53:03  mike
*	dev
*	
*	Revision 1.1  2001/05/05 23:53:51  mike
*	dev
*	
*
*/
package decodes.decoder;

import java.util.ArrayList;

import decodes.db.FormatStatement;

import ilex.var.TimedVariable;

/**
This is an abstract base class for all decodes operations.
*/
public abstract class DecodesOperation 
{
	/** List of timed variables that were generated by this operation */
	protected ArrayList<TimedVariable> decodedData = new ArrayList<TimedVariable>();
	
	/** Number of times to repeat this operation. */
	protected int repetitions;

	/** Position of the operation within the format statement */
	private TokenPosition tokenPosition = null;
	
	/** The format statement from which this operation was sourced */
	protected FormatStatement formatStatement = null;

	/**
	  Constructor is only called from subclass.
	  @param repetitions Number of times to repeat this operation
	*/
	protected DecodesOperation(int repetitions)
	{
		this.repetitions = repetitions;
	}

	/**
	  Each sub-class has a unique type-identifier.
	  @return unique type-identifier for concrete sub class.
	*/
	public abstract char getType();

	/**
	  @return Number of times to repeat this operation
	*/
	public int getRepetitions() { return repetitions; }

	/**
	  Subclass overides this operation to execute the operation against
	  the raw data contained in the passed DataOperations object, and
	  place the resulting data samples in the passed DecodedMessage object.

	  @param dd the DataOperations holds raw data and context.
	  @param msg the DecodedMessage into which decoded data is placed.

	  @throws DecoderException subclasses on various error conditions.
	*/
	public abstract void execute(DataOperations dd, DecodedMessage msg) 
		throws DecoderException;

	/**
	 * Deletes all stored results from previous executions.
	 */
	public void reset()
	{
		decodedData.clear();
	}
	
	/**
	 * Return the list of all timed variables resulting from previous execution.
	 * @return list of decoded timed variables.
	 */
	public ArrayList<TimedVariable> getDecodedData()
	{
		return decodedData;
	}

	public TokenPosition getTokenPosition()
	{
		return tokenPosition;
	}

	public void setTokenPosition(TokenPosition tokenPosition)
	{
		this.tokenPosition = tokenPosition;
	}

	public FormatStatement getFormatStatement()
	{
		return formatStatement;
	}

	public void setFormatStatement(FormatStatement formatStatement)
	{
		this.formatStatement = formatStatement;
	}
	
	

}

